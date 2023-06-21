/* scale.h */
/* header file for scale.c */

#define ANSI_SYS        1        /* 1: compile for ANSI_SYS color; 0: don't */

#define MAXL            256

#define K               30       /* points in the band-pass filter           */
#define HEAD_LINES      12       /* lines in the header of raw data file     */
#define BURST_RATE      1250.0   /* skew on HPADDA  is 750 micro-sec         */
#define DIGITAL_RANGE   0x7FFFFF /* 24 bit                                   */

//#define BURST_RATE    100000.0 /* skew on DAS1602 is  50 micro-sec         */
//efine DIGITAL_RANGE   0x8000   /* 16 bit                                   */

#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937510
#endif

void read_header ( char *argv[],
        int *nScan, float *sr, int *first_chnl, int *last_chnl, int *nChnl, float *range);

void read_sensitivity ( char *argv[], char *title,
		float *sensi, int *C, int *D, float *S,
		int start_chnl, int stop_chnl,
		char *Xlabel, char *Ylabel,
		char (*chnl_label)[MAXL], char (*units)[MAXL],
        	int *integ_chnl, int *diffr_chnl );

