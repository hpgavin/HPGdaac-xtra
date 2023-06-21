/* Program square.c - makes a file of an asymmetric square wave */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define REP_PNTS 1	/* repeated points (to reduce bandwidth) */

int main()
{

char	fn[80];		/* the output sinesweep filename	*/

FILE	*fp;		/* file pointer to the output data file	*/

float	t_pos, t_neg, t_zero, time, dt = 0.001, sr = 1000;

int	zero = 0,	/* output value for zero volts		*/
	pos  = 2000,	
	neg  = -100,
	pts=0, p, points,sv; /* number of points in the output file	*/

	/*	 Interactive  Input	*/

	printf(" Length of time of data: ");	sv=scanf("%f", &time);
	printf(" Samples per second: ");	sv=scanf("%f", &sr);
	printf(" Positive Amplitude: ");	sv=scanf("%d", &pos);
	printf(" Negative Amplitude: ");	sv=scanf("%d", &neg);
	printf(" Positive Time: ");		sv=scanf("%f", &t_pos);
	printf(" Negative Time: ");		sv=scanf("%f", &t_neg);
	printf(" Zero Time: ");			sv=scanf("%f", &t_zero);
	printf(" Output filename: ");		sv=scanf("%s", fn);

	dt = 1.0/sr;

	if (( fp = fopen(fn, "w" )) == NULL ) {
		printf ("  error: cannot open %s \n", fn );
		exit(1);
	}
	
	fprintf(fp, "%% positive amplitude = %d  ", pos );
	fprintf(fp, "%% positive time = %f \n", t_pos );
	fprintf(fp, "%% negative amplitude = %d  ", neg );
	fprintf(fp, "%% negative time = %f \n", t_neg );
	fprintf(fp, "%% zero time = %f \n", t_zero );
	fprintf(fp, "%% total time = %f  ", time );
	fprintf(fp, "%% time step = %f  \n", dt );
	fprintf(fp, "%5d\n%5d\n", zero, zero );
	pts += 3;

	while ( pts < (int)(time/dt) ) {
		for (p=1;p<(int)(t_zero/dt);p++) {
			fprintf(fp, "%5d\n", zero );
			pts++;
		}
		for (p=1;p<(int)(t_pos/dt);p++) {
			fprintf(fp, "%5d\n", pos);
			pts++;
		}
		for (p=1;p<(int)(t_neg/dt);p++) {
			fprintf(fp, "%5d\n", neg);
			pts++;
		}
	}
	fprintf( fp, "%5d\n%5d\n%5d\n", zero, zero, zero);
	pts += 3;
	printf("  %d data points in %f seconds\n", pts, (float)(pts*dt) );
	fclose (fp);
}
