/*************************************************************************
 Program skew.c  -  calculate the channel to channel skew of an AD device
 H.P. Gavin, Dept. of Civil Engineering, Univ. of Michigan, 28 May 93
**************************************************************************/

#include <stdio.h>	/* standard input/output library routines	*/
#include <stdlib.h>	/* standard input/output library routines	*/
#include <math.h>	/* math funcion library for fabs()		*/

#define	HEAD_LINES 12 

FILE	*fp;			/* pointer to the data file		*/

int main ( int argc, char *argv[] )
{
	char	ch = 'a';	/* indicates End Of File (EOF)		*/

	float	data, x=0, y=0, x_old=0,   /* data points for x and y	*/
		amp = 0, freq, chan_diff,
		dt, Sum_dt = 0, Sum_dt_dt = 0, dt_sd,
		sr, time=0, time_old=0,
		ZX = 0,		/* total number of zero crossings	*/
		sign(float x);	/* returns the sign of x */

	int	Xcol, Ycol,	/* columns for X and Y data		*/
		Ncol,		/* number of columns in the data file	*/
		points = 0,	/* total number of data points		*/
		i;

	void	usr_inp();	/* open input file, read control data	*/

	usr_inp ( argv, &Ncol, &Xcol, &Ycol, &sr );

	chan_diff = fabs( (float) Xcol - (float) Ycol );

	while ( ch != EOF ) {
		for (i=1; i<= Ncol; i++) {
			ch=fscanf(fp,"%f", &data);
			if (ch == EOF) break;
			if (i==Xcol) x=data;
			if (i==Ycol) y=data;
		}
		if (ch == EOF) break;
		dt = ( y - x ) * (x-x_old) / fabs(x-x_old);
		Sum_dt += dt;
		Sum_dt_dt += dt*dt;
		if (fabs(y) > amp) amp = fabs(y);
		if (fabs(x) > amp) amp = fabs(x);
		if ( x_old < 0 && x_old*x < 0 ) {
			if (ZX) time += points/sr - time_old;
			time_old = points/sr;
			ZX++;
		}
		++points;
		printf("%4d  %12.1f   %12.1f   %3.1f\n",
			points, dt, (x-x_old),  (x-x_old)/fabs(x-x_old) );
		x_old = x;
	}

	freq = ZX / time;

	dt_sd = sqrt( (Sum_dt_dt - Sum_dt*Sum_dt/points) / (points - 1.0) );

	fprintf(stderr,"%4d  points\n", points);
	fprintf(stderr,"  skew= %4.1f micro-sec  C.o.V.= %f",
	1e6*Sum_dt/(points*4.*amp*freq*chan_diff), dt_sd*points/Sum_dt );
	fprintf(stderr,"  max. amplitude= %.1f\n", amp );
	fprintf(stderr,"  triangle wave frequency= %.2f Hz\n", freq );
}

/*---------------------------------------------------------------------------
PROCEDURES	USR_INP  -  open input file, read in control data  25jun93
----------------------------------------------------------------------------*/

void usr_inp ( char *argv[], int *Ncol, int *Xcol, int *Ycol, float *sr )
{
	char	ch;
	int	i, sv;

	if ((fp = fopen(argv[1],"r")) == NULL ) {
		fprintf(stderr,"  error: cannot open file '%s'\n", argv[1]);
		fprintf(stderr,"  usage: skew [input file]\n");
		exit(0);
	}

	fprintf(stderr," Number of Columns: ");		sv=scanf("%d", Ncol);
	fprintf(stderr," X-data column: ");		sv=scanf("%d", Xcol);
	fprintf(stderr," Y-data column: ");		sv=scanf("%d", Ycol);
	fprintf(stderr," Sample Rate (Hz): ");		sv=scanf("%f", sr );

	for (i=1; i<= HEAD_LINES; i++) while((ch=fgetc(fp)) != '\n') ;

	return;
}

float sign ( float x )
{
	float	answer = 1.0;
	if ( x < 0.0 )	answer = -1.0;

	printf(" ... %f   %f ...", x, answer );

	return answer;
}
