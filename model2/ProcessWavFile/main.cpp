
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "WAVheader.h"
#include "common.h"
#include "stdfix_emu.h"
#include "fixed_point_math.h"



DSPfract x_history0[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };	//LEFT SURROUND
DSPfract y_history0[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };

DSPfract x_history1[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };	//LEFT
DSPfract y_history1[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };

DSPfract x_history2[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };	//BASS
DSPfract y_history2[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };

DSPfract x_history3[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };	//RIGHT SURROUND
DSPfract y_history3[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };

DSPfract x_history4[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };	//RIGHT
DSPfract y_history4[2] = { FRACT_NUM(0.0), FRACT_NUM(0.0) };

DSPfract filter2low[6] = { FRACT_NUM(0.07752551285593004), FRACT_NUM(0.15505102571186008), FRACT_NUM(0.07752551285593004), FRACT_NUM(0.50000000000000000000), FRACT_NUM(-0.310102051443364335), FRACT_NUM(0.120204102886728765)};
DSPfract filter2high[6] = { FRACT_NUM(0.37855335859929762), FRACT_NUM (-0.75710671719859525), FRACT_NUM(0.37855335859929762), FRACT_NUM(0.50000000000000000000), FRACT_NUM (-0.72712179314438065), FRACT_NUM(0.28703095755542768)};

DSPfract sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];
DSPfract tempBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];
DSPfract tempLeft[BLOCK_SIZE];
DSPfract tempRight[BLOCK_SIZE];
char decibels[50];
char* pEnd;

DSPfract left_s_channel[BLOCK_SIZE];
DSPfract left_channel[BLOCK_SIZE];
DSPfract bass_channel[BLOCK_SIZE];
DSPfract right_s_channel[BLOCK_SIZE];
DSPfract right_channel[BLOCK_SIZE];

DSPfract tempLeftHelp[BLOCK_SIZE];
DSPfract tempRightHelp1[BLOCK_SIZE];
DSPfract tempRightHelp2[BLOCK_SIZE];

enum output_mode { MODE_2_0_0, MODE_2_2_1 };
output_mode outputMode = MODE_2_0_0;
DSPint enable = 1; //1=on, 0=off
DSPfract input_gain = FRACT_NUM(1.0);

DSPfract gain_9 = FRACT_NUM(0.33884415613920255);
DSPfract gain_3 = FRACT_NUM(0.6456542290346555);

DSPaccum second_order_IIR(DSPfract input, DSPfract* coefficients, DSPfract* x_history, DSPfract* y_history);
void processing();
//float dBToinput_gain();

DSPint main(int argc, char* argv[])
{
	FILE *wav_in = NULL;
	FILE *wav_out = NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr, outputWAVhdr;

	// Init channel buffers
	/*for (DSPint i = 0; i < MAX_NUM_CHANNEL; i++)
		//memset(&sampleBuffer[i], 0, BLOCK_SIZE);
		sampleBuffer[i][j] = FRACT_NUM(0.0);*/
	for (DSPint i = 0; i < MAX_NUM_CHANNEL; i++)
	{
		for (DSPint j = 0; j < BLOCK_SIZE; j++) {
			sampleBuffer[i][j] = FRACT_NUM(0.0);
			tempBuffer[i][j] = FRACT_NUM(0.0);
		}
	}
	if (argc != 6)
	{
		printf("Nema dovoljno argumenata!\n");
		return -1;
	}

	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName, argv[1]);
	wav_in = OpenWavFileForRead(WavInputName, "rb");
	strcpy(WavOutputName, argv[2]);
	wav_out = OpenWavFileForRead(WavOutputName, "wb");

	//mode selection 
	enable = atoi(argv[3]);
	printf("Enable: 1-ON, 0-OFF -> %d\n", enable);

	//strcpy(decibels, argv[4]);
	//input_gain = atof(decibels);
	//input_gain = dBToinput_gain();
	printf("Gain: %f dB\n", input_gain.toLong());

	input_gain = FRACT_NUM(atof(argv[4]));

	int outputMode1 = atoi(argv[5] + 4);
	if (outputMode1 == 0)
	{
		printf("Mode: 2_0_0!\n");
		outputMode = MODE_2_0_0;
	}
	else
	{
		printf("Mode: Default 2_2_1\n");
		outputMode = MODE_2_2_1;
	}

	//-------------------------------------------------

	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in, inputWAVhdr);
	//-------------------------------------------------

	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	//outputWAVhdr.fmt.NumChannels = inputWAVhdr.fmt.NumChannels; // change number of channels
	outputWAVhdr.fmt.NumChannels = NUM_OUT_CHANNELS;
	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size / inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate / inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign / inputWAVhdr.fmt.NumChannels;

	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out, outputWAVhdr);


	// Processing loop
	//-------------------------------------------------	
	{
		DSPint sample;
		DSPint BytesPerSample = inputWAVhdr.fmt.BitsPerSample / 8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		DSPint iNumSamples = inputWAVhdr.data.SubChunk2Size / (inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample / 8);

		// exact file length should be handled correctly...
		for (DSPint i = 0; i < iNumSamples / BLOCK_SIZE; i++)
		{
			for (DSPint j = 0; j < BLOCK_SIZE; j++)
			{
				for (DSPint k = 0; k < inputWAVhdr.fmt.NumChannels; k++)
				{
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);

					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;
					tempBuffer[k][j] = sample / SAMPLE_SCALE;	// scale sample to 1.0/-1.0 range		

				}
			}


			processing();


			for (DSPint j = 0; j<BLOCK_SIZE; j++)
			{
				for (DSPint k = 0; k<outputWAVhdr.fmt.NumChannels; k++)
				{
					sample = sampleBuffer[k][j].toLong();	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample / 8, 1, wav_out);
				}
			}
		}
	}

	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}

DSPaccum second_order_IIR(DSPfract input, DSPfract* coefficients, DSPfract* x_history, DSPfract* y_history) {
	DSPaccum output = 0;


	output += *coefficients * input << 1;        /* A0 * x(n)     */
	output += (*(coefficients + 1) * *x_history) << 1; /* A1 * x(n-1) */
	output += (*(coefficients + 2) * *(x_history + 1)) << 1; /* A2 * x(n-2)   */
	output -= (*(coefficients + 4) * *y_history) << 1; /* B1 * y(n-1) */
	output -= (*(coefficients + 5) * *(y_history + 1)) << 1; /* B2 * y(n-2)   */

	*(y_history + 1) = *y_history;    /* y(n-2) = y(n-1) */
	*y_history = output; /* y(n-1) = y(n)   */
	*(x_history + 1) = *x_history;  /* x(n-2) = x(n-1) */
	*x_history = input;          /* x(n-1) = x(n)   */
	


	return output;
}
void processing()
{
	DSPint i;
	DSPint k;
	DSPfract tempLeft, tempRight;

	DSPfract* left_s_channel = sampleBuffer[0];
	DSPfract* left_channel = sampleBuffer[1];
	DSPfract* bass_channel = sampleBuffer[2];
	DSPfract* right_s_channel = sampleBuffer[3];
	DSPfract* right_channel = sampleBuffer[4];

	DSPfract* tempLeftHelp = tempBuffer[0];
	DSPfract* tempRightHelp1 = tempBuffer[1];
	DSPfract* tempRightHelp2 = tempBuffer[2];

	

	for (i = 0; i < BLOCK_SIZE; i++)
	{
		tempLeft = (*left_s_channel) * input_gain;
		tempRight = (*left_channel) * input_gain;
		
		if (enable == 1)
		{
			switch (outputMode)
			{
			case MODE_2_2_1:
				/*LEFT CHANNEL*/
				for (k = 0; k < 2; k++)
				{
					*left_s_channel = second_order_IIR(tempLeft, filter2high, x_history0, y_history0);
				}
				for (k = 0; k < 2; k++)
				{
					*bass_channel = second_order_IIR(tempLeft, filter2low, x_history2, y_history2);
				}
				*bass_channel = (*bass_channel) * gain_9;
				for (k = 0; k < 2; k++)
				{
					*tempLeftHelp = second_order_IIR(*left_s_channel, filter2low, x_history1, y_history1);
					
				}
				*tempLeftHelp = (*tempLeftHelp) * gain_3;
				*left_channel = (*tempLeftHelp).toDouble() + (*bass_channel).toDouble();
				/*RIGHT CHANNEL*/

				for (k = 0; k < 2; k++)
				{
					*right_s_channel = second_order_IIR(tempRight, filter2high, x_history3, y_history3);
				}
				for (k = 0; k < 2; k++)
				{
					*tempRightHelp1 = second_order_IIR(*right_s_channel, filter2low, x_history4, y_history4);
				}
				*tempRightHelp1 = (*tempRightHelp1) * gain_3;
				for (k = 0; k < 2; k++)
				{
					*tempRightHelp2 = second_order_IIR(tempRight, filter2low, x_history4, y_history4);
				}
				*tempRightHelp2 = (*tempRightHelp2) * gain_9;
				*right_channel = (*tempRightHelp1) + (*tempRightHelp2);
				
				break;

			case MODE_2_0_0:
				/*LEFT CHANNEL*/
				for (k = 0; k < 2; k++)
				{
					*left_s_channel = second_order_IIR(tempLeft, filter2high, x_history0, y_history0);
				}
				for (k = 0; k < 2; k++)
				{
					*bass_channel = second_order_IIR(tempLeft, filter2low, x_history2, y_history2);
				}
				*bass_channel = *bass_channel * gain_9;
				for (k = 0; k < 2; k++)
				{
					*tempLeftHelp = second_order_IIR(*left_s_channel, filter2low, x_history1, y_history1);
				}
				*tempLeftHelp = *tempLeftHelp * gain_3;
				*left_channel = *tempLeftHelp + *bass_channel;
				*left_s_channel = *left_s_channel * FRACT_NUM(0.0);
				*bass_channel = *bass_channel * FRACT_NUM(0.0);

				/*RIGHT CHANNEL*/

				for (k = 0; k < 2; k++)
				{
					*right_s_channel = second_order_IIR(tempRight, filter2high, x_history3, y_history3);
				}
				for (k = 0; k < 2; k++)
				{
					*tempRightHelp1 = second_order_IIR(*right_s_channel, filter2low, x_history4, y_history4);
				}
				*tempRightHelp1 = *tempRightHelp1 * gain_3;
				for (k = 0; k < 2; k++)
				{
					*tempRightHelp2 = second_order_IIR(tempRight, filter2low, x_history4, y_history4);
				}
				*tempRightHelp2 = *tempRightHelp2 * gain_9;
				*right_channel = *tempRightHelp1 + *tempRightHelp2;
				*right_s_channel = *right_s_channel * FRACT_NUM(0.0);
			
			}
		}
		else
		{
			/*LEFT CHANNEL*/
			 
			*left_s_channel = FRACT_NUM(0.0);
			*bass_channel = FRACT_NUM(0.0);
			*left_channel = FRACT_NUM(0.0);

			/*RIGHT CHANNEL*/

			*right_s_channel = FRACT_NUM(0.0);			
			*right_channel = FRACT_NUM(0.0);
		}
		left_s_channel++;
		left_channel++;
		bass_channel++;
		right_s_channel++;
		right_channel++;
		tempLeftHelp++;
		tempRightHelp1++;
		tempRightHelp2++;

		
	}
}
/*float dBToinput_gain()
{
	return pow(10.0f, input_gain/20.0f);
}*/