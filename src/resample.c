/*******************************************************************************
resample.c -  Re-samples a multi-column time history data file.
	To compile:	gcc -O -o resample resample.c  -lm
	To run:		resample in_file  out_file
H.P. Gavin, Civil and Environmental Engineering, Duke University, 30 March 1999 
*******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "../../HPGnumlib/HPGutil.h"

int main ( int argc, char *argv[] )
{
	char	cha='a', chb='a', line[MAXL], integers;

	FILE	*fp_in,		/*  input file pointer			*/
		*fp_out;	/* output file pointer			*/

	float	x[20], x_old[20], max[20], xi, i; 
	int	points=0, ok=1, head_lines = 0, c, j, len, N, Ncol, sv;


	if ((fp_in = fopen (argv[1], "r")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[1]);
		printf (" usage: resample infile outfile\n");
		exit(0);
	}
	if ((fp_out = fopen (argv[2], "w")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[1]);
		printf (" usage: resample infile outfile\n");
		exit(0);
	}

	printf("  Resampleing factor (integer>1): ");	sv=scanf("%d", &N); 
	printf("  Number of Columns: ");		sv=scanf("%d", &Ncol);
	(void) getchar();
	printf("  Output integers? (y/n): ");		sv=scanf("%c", &integers); 

				/* determine the size of the header */
	do {
		cha = getc(fp_in);
		if ( cha == '#' || cha == '%' ) {
			++head_lines;
			while (( chb = getc(fp_in) ) != '\n' ) ;
		}
	} while ( cha == '#' || cha == '%' );
	rewind (fp_in);

	for (c=1; c<=Ncol; c++)	max[c] = 0;
	for (c=1; c<=Ncol; c++)	x_old[c] = 0.0;

	for (i=1; i<=head_lines; i++) 		/* scan through the header */
		len = getLine ( fp_in, MAXL, line );

	while ( ok == 1 ) {			/* find maximum values	*/
		for (c=1; c<=Ncol; c++)	 ok = fscanf(fp_in,"%f", &x[c]);
		for (c=1; c<=Ncol; c++)	
			if ( fabs(x[c]) > max[c] )   max[c] = fabs(x[c]);
	}
	rewind (fp_in);	

	for (i=1; i<=head_lines; i++) {		/* copy the header */
		len = getLine ( fp_in, MAXL, line );
		for (j=0; j<len; j++)	putc ( line[j], fp_out );
	}

	ok = 1;
	while ( ok == 1 ) {

		for (c=1; c<=Ncol; c++)		ok = fscanf(fp_in,"%f", &x[c]);

		for (i=0; i<N; i++) {
		    for (c=1; c<=Ncol; c++) { /* linear interpolation */
			xi = x_old[c] + (x[c]-x_old[c])*(float)i/N;
			if (integers == 'y')
				fprintf(fp_out,"%6.0f", xi*2000/max[c] );
			else
				fprintf(fp_out,"%15.6e", xi );
		    }
		    fprintf(fp_out,"\n");
		}
		for (c=1; c<=Ncol; c++)		x_old[c] = x[c];

		++points;
	}

	printf("   %d points\n", points );

	fclose (fp_in);
	fclose (fp_out);
}


