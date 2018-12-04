
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "WAVheader.h"
#include "iir_filter.h"

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8
#define NUM_OUT_CHANNELS 5
#define GAIN_4DB 0.6309573444801932
#define GAIN_9_4 0.33884415613920255
#define GAIN_3_8 0.6456542290346555

double x_history0[2] = { 0.0, 0.0 };	//LEFT SURROUND
double y_history0[2] = { 0.0, 0.0 };

double x_history1[2] = { 0.0, 0.0 };	//LEFT
double y_history1[2] = { 0.0, 0.0 };

double x_history2[2] = { 0.0, 0.0 };	//BASS
double y_history2[2] = { 0.0, 0.0 };

double x_history3[2] = { 0.0, 0.0 };	//RIGHT SURROUND
double y_history3[2] = { 0.0, 0.0 };

double x_history4[2] = { 0.0, 0.0 };	//RIGHT
double y_history4[2] = { 0.0, 0.0 };


double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];
double tempBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];
char decibels[50];
char* pEnd;

enum output_mode  {MODE_2_0_0, MODE_2_2_1}; 
output_mode outputMode = MODE_2_0_0;
int enable = 1; //1=on, 0=off
float input_gain = 1.0;

double second_order_IIR(double input, double* coefficients, double* x_history, double* y_history);
void processing();
float dBToinput_gain();

int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

	// Init channel buffers
	for(int i=0; i<MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i],0,BLOCK_SIZE);

	if (argc != 6)
	{
		printf("Nema dovoljno argumenata!\n");
		return -1;
	}

	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName,argv[1]);
	wav_in = OpenWavFileForRead (WavInputName,"rb");
	strcpy(WavOutputName,argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName,"wb");

	//mode selection 
	enable = atoi(argv[3]);
	printf("Enable: 1-ON, 0-OFF -> %d\n", enable);

	strcpy(decibels, argv[4]);
	input_gain = strtol(decibels, &pEnd, 10);
	input_gain = dBToinput_gain(); 
	printf("Gain: %f dB\n", input_gain);

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
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	//outputWAVhdr.fmt.NumChannels = inputWAVhdr.fmt.NumChannels; // change number of channels
	outputWAVhdr.fmt.NumChannels = NUM_OUT_CHANNELS;
	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);


	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
		// exact file length should be handled correctly...
		for (int i = 0; i < iNumSamples / BLOCK_SIZE; i++)
		{
			for (int j = 0; j < BLOCK_SIZE; j++)
			{
				for (int k = 0; k < inputWAVhdr.fmt.NumChannels; k++)
				{
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);

					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;	
					//tempBuffer[k][j] = sample / SAMPLE_SCALE;	// scale sample to 1.0/-1.0 range		

				}
			}


			processing();
		

			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j] * SAMPLE_SCALE ;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample/8, 1, wav_out);		
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

double second_order_IIR(double input, double* coefficients, double* x_history, double* y_history) {
	double output = 0;

	output += coefficients[0] * input;        /* A0 * x(n)     */
	output += coefficients[1] * x_history[0]; /* A1 * x(n-1) */
	output += coefficients[2] * x_history[1]; /* A2 * x(n-2)   */
	output -= coefficients[4] * y_history[0]; /* B1 * y(n-1) */
	output -= coefficients[5] * y_history[1]; /* B2 * y(n-2)   */


	y_history[1] = y_history[0];    /* y(n-2) = y(n-1) */
	y_history[0] = output; /* y(n-1) = y(n)   */
	x_history[1] = x_history[0];  /* x(n-2) = x(n-1) */
	x_history[0] = input;          /* x(n-1) = x(n)   */

	return output;
}
void processing()
{
	int i;
	int k;

	double tempL, tempR;

	for (i = 0; i < BLOCK_SIZE; i++)
	{
		tempL = sampleBuffer[0][i] * input_gain;
		tempR = sampleBuffer[1][i] * input_gain;

		switch(outputMode)
		{
		case MODE_2_0_0:
			if (enable == 1) {
				/*LEFT CHANNEL*/
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[0][i] = second_order_IIR(tempL, filter2high, x_history0, y_history0);
				}
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[2][i] = second_order_IIR(tempL, filter2low, x_history2, y_history2);
				}
				sampleBuffer[2][i] = sampleBuffer[2][i] * GAIN_9_4;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[0][i] = second_order_IIR(sampleBuffer[0][i], filter2low, x_history1, y_history1);
				}
				tempBuffer[0][i] = tempBuffer[0][i] * GAIN_3_8;
				sampleBuffer[1][i] = tempBuffer[0][i] + sampleBuffer[2][i];
				sampleBuffer[0][i] = sampleBuffer[0][i] * 0.0;
				sampleBuffer[2][i] = sampleBuffer[2][i] * 0.0;

				/*RIGHT CHANNEL*/

				for (k = 0; k < 2; k++)
				{
					sampleBuffer[3][i] = second_order_IIR(tempR, filter2high, x_history3, y_history3);
				}
				for (k = 0; k < 2; k++)
				{
					tempBuffer[1][i] = second_order_IIR(sampleBuffer[3][i], filter2low, x_history4, y_history4);
				}
				tempBuffer[1][i] = tempBuffer[1][i] * GAIN_3_8;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[2][i] = second_order_IIR(tempR, filter2low, x_history4, y_history4);
				}
				tempBuffer[2][i] = tempBuffer[2][i] * GAIN_9_4;
				sampleBuffer[4][i] = tempBuffer[1][i] + tempBuffer[2][i];
				sampleBuffer[3][i] = sampleBuffer[3][i] * 0.0;
			}
			else 
			{
				/*LEFT CHANNEL*/
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[0][i] = second_order_IIR(tempL, filter2high, x_history0, y_history0) * 0.0;
				}
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[2][i] = second_order_IIR(tempL, filter2low, x_history2, y_history2);
				}
				sampleBuffer[2][i] = sampleBuffer[2][i] * GAIN_9_4 * 0.0;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[0][i] = second_order_IIR(sampleBuffer[0][i], filter2low, x_history1, y_history1);
				}
				tempBuffer[0][i] = tempBuffer[0][i] * GAIN_3_8;
				sampleBuffer[1][i] = tempBuffer[0][i] + sampleBuffer[2][i] * 0.0;

				/*RIGHT CHANNEL*/

				for (k = 0; k < 2; k++)
				{
					sampleBuffer[3][i] = second_order_IIR(tempR, filter2high, x_history3, y_history3) * 0.0;
				}
				for (k = 0; k < 2; k++)
				{
					tempBuffer[1][i] = second_order_IIR(sampleBuffer[3][i], filter2low, x_history4, y_history4);
				}
				tempBuffer[1][i] = tempBuffer[1][i] * GAIN_3_8;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[2][i] = second_order_IIR(tempR, filter2low, x_history4, y_history4);
				}
				tempBuffer[2][i] = tempBuffer[2][i] * GAIN_9_4;
				sampleBuffer[4][i] = tempBuffer[1][i] + tempBuffer[2][i] * 0.0;
			}
			break;

		case MODE_2_2_1:
			if (enable == 1) {
				/*LEFT CHANNEL*/
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[0][i] = second_order_IIR(tempL, filter2high, x_history0, y_history0);
				}
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[2][i] = second_order_IIR(tempL, filter2low, x_history2, y_history2);
				}
				sampleBuffer[2][i] = sampleBuffer[2][i] * GAIN_9_4;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[0][i] = second_order_IIR(sampleBuffer[0][i], filter2low, x_history1, y_history1);
				}
				tempBuffer[0][i] = tempBuffer[0][i] * GAIN_3_8;
				sampleBuffer[1][i] = tempBuffer[0][i] + sampleBuffer[2][i];

				/*RIGHT CHANNEL*/

				for (k = 0; k < 2; k++)
				{
					sampleBuffer[3][i] = second_order_IIR(tempR, filter2high, x_history3, y_history3);
				}
				for (k = 0; k < 2; k++)
				{
					tempBuffer[1][i] = second_order_IIR(sampleBuffer[3][i], filter2low, x_history4, y_history4);
				}
				tempBuffer[1][i] = tempBuffer[1][i] * GAIN_3_8;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[2][i] = second_order_IIR(tempR, filter2low, x_history4, y_history4);
				}
				tempBuffer[2][i] = tempBuffer[2][i] * GAIN_9_4;
				sampleBuffer[4][i] = tempBuffer[1][i] + tempBuffer[2][i];
			}
			else
			{
				/*LEFT CHANNEL*/
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[0][i] = second_order_IIR(tempL, filter2high, x_history0, y_history0) * 0.0;
				}
				for (k = 0; k < 2; k++)
				{
					sampleBuffer[2][i] = second_order_IIR(tempL, filter2low, x_history2, y_history2);
				}
				sampleBuffer[2][i] = sampleBuffer[2][i] * GAIN_9_4 * 0.0;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[0][i] = second_order_IIR(sampleBuffer[0][i], filter2low, x_history1, y_history1);
				}
				tempBuffer[0][i] = tempBuffer[0][i] * GAIN_3_8;
				sampleBuffer[1][i] = tempBuffer[0][i] + sampleBuffer[2][i] * 0.0;

				/*RIGHT CHANNEL*/

				for (k = 0; k < 2; k++)
				{
					sampleBuffer[3][i] = second_order_IIR(tempR, filter2high, x_history3, y_history3) * 0.0;
				}
				for (k = 0; k < 2; k++)
				{
					tempBuffer[1][i] = second_order_IIR(sampleBuffer[3][i], filter2low, x_history4, y_history4);
				}
				tempBuffer[1][i] = tempBuffer[1][i] * GAIN_3_8;
				for (k = 0; k < 2; k++)
				{
					tempBuffer[2][i] = second_order_IIR(tempR, filter2low, x_history4, y_history4);
				}
				tempBuffer[2][i] = tempBuffer[2][i] * GAIN_9_4;
				sampleBuffer[4][i] = tempBuffer[1][i] + tempBuffer[2][i] * 0.0;
			}
		}

	}
}
float dBToinput_gain()
{
	return pow(10.0f, input_gain / 20.0f);
}