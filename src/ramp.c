/******************************************************************************
ramp.c - makes a file of ramp data with constant velocity ... 

Henri Gavin, Dept. Civil Engineering, Duke University, 22 Dec. 2003
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

int main()
{
	char	fn[80],	// the output data filename
		ch;

	FILE	*f;	/* file pointer to the output data file	*/

	float	sr, n,		/* sample rate		*/
		zero = 0,	/* output value for zero volts		*/
		ampl = 2000,	/* amplitude				*/
		ramp_time,	/* time for the ramp	*/
		pre_ramp,	/* time for the ramp	*/
		post_ramp;	/* time for the ramp	*/

	int	points = 0, sv;


	/*	 Interactive  Input	*/

	printf(" sample rate: ");		sv=scanf("%f", &sr);
	printf(" ramp time (s): ");		sv=scanf("%f", &ramp_time);
	printf(" pre ramp time (s): ");		sv=scanf("%f", &pre_ramp);
	printf(" post ramp time (s): ");	sv=scanf("%f", &post_ramp);

	printf(" Output ramp waveform filename: ");	sv=scanf("%s", fn);

	if (( f = fopen(fn, "w" )) == NULL ) {
		printf ("  error: cannot open %s \n", fn );
		exit(0);
	}
	
	fprintf(f, "%5.0f\n%5.0f\n", zero, zero );

	for (n=0; n < pre_ramp*sr; n++) {
		fprintf(f, "%5.0f\n", zero );
		points++;
	}
	for (n=0; n < ramp_time*sr; n++) {
		fprintf(f, "%5.0f\n", ampl*(zero+n/sr/ramp_time) );
		points++;
	}
	for (n=0; n < post_ramp*sr; n++) {
		fprintf(f, "%5.0f\n", ampl );
		points++;
	}
	for (n=0; n < 100*ramp_time*sr; n++) {
		fprintf(f, "%5.0f\n", ampl*(1-(zero+n/sr/(100*ramp_time)) ));
		points++;
	}
	fprintf(f, "%5.0f\n%5.0f", zero, zero );
	printf("  Test time = %f \n", points/sr ); 

	fclose(f);
}
