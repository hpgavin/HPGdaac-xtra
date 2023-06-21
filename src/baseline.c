/*******************************************************************************
baseline.c - Perform baseline correction on waveform data
		to complile:	gcc -O -o baseline baseline.c
		to run:		baseline infile outfile
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../HPGnumlib/HPGutil.h"

int main ( int argc, char *argv[])
{
	FILE	*fp_in,		//  input file pointer
		*fp_out;	// output file pointer

        char    c = 'a', ch = 'a', line[MAXL]; 

	float	x[20], d_o, d_f, SF, Po, Pf, peak=0.0;
	int	i, j, Ncol, Col, ok = 1,
		head_lines=0, points=0, p=0, len, sfrv = 0,
                getLine ( FILE *fp, int lim, char *s );


	if ((fp_in = fopen (argv[1], "r")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[1]);
		printf (" usage: baseline infile outfile\n");
		exit(0);
	}
	if ((fp_out = fopen (argv[2], "w")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[2]);
		printf (" usage: baseline infile outfile\n");
		exit(0);
	}

	printf(" Number of columns: ");		sfrv = scanf("%d", &Ncol );
	printf(" Column to correct: ");		sfrv = scanf("%d", &Col ); 
	printf(" Scaling factor: ");	  	sfrv = scanf("%f", &SF ); 
	printf(" Number of start-points to average: ");	sfrv = scanf("%f", &Po );
	printf(" Number of end-points to average: ");	sfrv = scanf("%f", &Pf );
	if (Po <= 0)  Po = 1.0;
	if (Pf <= 0)  Pf = 1.0;

			/* find the size of the header	*/
	do {
		ch = getc(fp_in);
		if ( ch == '#' || ch == '%' ) {
			++head_lines;
			while (( c = getc(fp_in) ) != '\n' ) ;
		}
	} while ( ch == '#' || ch == '%' );
	rewind (fp_in);

	d_o = 0.0;
	for (i=1;i<=head_lines;i++)     /* scan through the header */
		while (( ch = getc(fp_in)) != '\n') ;

	while ( ok == 1 ) {
		++points; 
		for ( i = 1; i <= Ncol; i++ )	
			ok = fscanf(fp_in,"%f", &x[i] );
		if ( points <= Po && ok == 1 )	d_o += x[Col];
		if ( peak < fabs(x[Col]) ) peak = fabs(x[Col]);
	}
	d_o /= Po;
	--points;
	
	printf(" %d points\n", points );


	rewind(fp_in);
	for (i=1;i<=head_lines;i++)     /* scan through the header */
		while (( ch = getc(fp_in)) != '\n') ;
	ok = 1;
	p = 0;
	d_f = 0.0;
	while ( ok == 1 ) { 
		++p; 
		for ( i = 1; i <= Ncol; i++ )	
			ok = fscanf(fp_in,"%f", &x[i] );  
		if ( p > points-Pf && ok == 1 )	d_f += x[Col];
	}  
	d_f /= Pf;

	rewind(fp_in);
        for (i=1; i<=head_lines; i++) {         /* copy the header */
                len = getLine ( fp_in, MAXL, line );
                for (j=0; j<len; j++) putc ( line[j], fp_out );
        }
	ok = 1;
	p  = 0;

	SF = 2040.0 / peak;
/*	SF = 1.0;		*/

	while ( ok == 1 ) {
		++p; 
		for ( i = 1; i <= Ncol; i++ )	
			ok = fscanf(fp_in,"%f", &x[i] );
		if ( ok == 1 ) 
			fprintf(fp_out," %7.0f \n", 
			  ( x[Col] - d_o - (d_f-d_o)*(p-1)/(points-1) ) * SF );
	}

	fclose (fp_in);
	fclose (fp_out);
}

