#ifndef COMMON_H
#define COMMON_H
#include "stdfix_emu.h"
#include "fixed_point_math.h"

/* Basic constants */

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8
#define NUM_OUT_CHANNELS 5
#define GAIN_4DB 0.6309573444801932
#define GAIN_9_4 0.33884415613920255
#define GAIN_3_8 0.6456542290346555

/*enum output_mode { MODE_2_0_0, MODE_2_2_1 };
output_mode outputMode = MODE_2_0_0;
DSPint enable = 1; //1=on, 0=off
DSPfract input_gain = FRACT_NUM(1.0);
char decibels[50];
char* pEnd;*/


/* DSP type definitions */
typedef short DSPshort;					/* DSP integer */
typedef unsigned short DSPushort;		/* DSP unsigned integer */
typedef int DSPint;						/* native integer */
typedef fract DSPfract;					/* DSP fixed-point fractional */
typedef long_accum DSPaccum;			/* DSP fixed-point fractional */

#endif
