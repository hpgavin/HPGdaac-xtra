/**************************************************************************
 Program limits.c  -  scans a data file and keeps track of max and min 
 H.P. Gavin, Dept. of Civil Engineering, Univ. of Michigan, 12 Mar 95
***************************************************************************/

#include <stdio.h>	/* standard input/output library routines	*/
#include <math.h>	/* standard input/output library routines	*/
#include <stdlib.h>

#define MAXL        128      /* longest header line length          */

FILE	*fp_in;			/* pointer to the data file		*/

int main ( int	argc, char *argv[] )
{
	int 	head_lines = 0,	/* number of lines in the header	*/
		ch = 1,		/* number of successfully read values	*/
		col_start, col_stop,	/* starting and stopping channels */
		points = 0,	/* total number of data points		*/
		col;

	float	data,			/* a data value			*/
		delta_t,		/* time step			*/
		Vs, Ds,			/* veloc. and displ. scale factors */
		max[16], min[16],	/* max and min values		*/
		Tmax[16], Tmin[16],	/* Times of max and min values	*/
		avg[16], rms[16];	/* average and rms values	*/

	void	usr_inp();	/* open input file, read control data	*/

	usr_inp ( argv, &head_lines, &col_start, &col_stop, &delta_t );

	for (col = 0; col <= 15; col++)	
	  avg[col] = rms[col] = max[col] = min[col] = Tmax[col] = Tmin[col] = 0;

	while ( ch == 1 ) {
		for (col = col_start; col <= col_stop; col++) {
			ch=fscanf(fp_in,"%f", &data);
			if (ch == EOF) break;
			if ( points == 0 ) max[col] = min[col] = data;
			avg[col] += data;
			rms[col] += data * data;
			if (data > max[col]) {
				max[col] = data; 
				Tmax[col] = points*delta_t;
			}
			if (data < min[col]) {
				min[col] = data; 
				Tmin[col] = points*delta_t;
			}
		}
		if (ch != 1) break;
		++points;
	}
	for ( col = col_start; col <= col_stop; col++ ) {
		avg[col] /= points;
		rms[col] = sqrt ( rms[col] / (float) points );
	}

	Vs =  (-min[2] > max[2]) ? -min[2] : max[2];
	Ds =  (-min[3] > max[3]) ? -min[3] : max[3];

	Vs = Vs / 50.0;
	Ds = Ds / 7.0;

	printf (" %s : %4d  points   Vs = %6.2f   Ds = %6.2f  ",
						argv[1], points, Vs, Ds);
	if (Vs < Ds)	printf(" displ. limited\n");
	else		printf(" veloc. limited\n");
	printf ("# ____________________________________________________________________\n");
	printf ("# column   minimum   maximum   average    RMS       Tmin       Tmax\n");

	for ( col = col_start; col <= col_stop; col ++)
		printf("# %3d     %8.3f  %8.3f  %8.3f  %8.3f   %8.3f   %8.3f\n",
			col, min[col], max[col], avg[col], rms[col], Tmin[col], Tmax[col] );
	printf ("# ____________________________________________________________________\n");
}

/*-----------------------------------------------------------------------------
USR_INP  -  open input file, read in control data 			12mar95
------------------------------------------------------------------------------*/

void usr_inp ( argv, head_lines, col_start, col_stop, delta_t )
char	*argv[];
float	*delta_t;
int	*col_start, *col_stop, *head_lines;
{
	char	ch, c='a';
	int	i;

	if ((fp_in = fopen(argv[1],"r")) == NULL ) {
		fprintf(stderr,"  error: cannot open file '%s'\n", argv[1]);
		fprintf(stderr,"  usage: skew [input file]\n");
		exit(1);
	}

	*col_start = 1;
	*col_stop  = 3;
	*delta_t   = 0.020;

                        /* find the size of the header  */
        do {
                if (( ch = getc(fp_in) ) == '#' ) {
                        ++(*head_lines);
                        while (( c = getc(fp_in) ) != '\n' ) ;
                }
        } while ( ch == '#' );
        rewind (fp_in);

        for (i=1;i<= *head_lines; i++)     /* scan through the header */
                while (( ch = getc(fp_in)) != '\n') ;

	return;
}
