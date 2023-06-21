/************************************************************************/
/* Program deskip.c - resorts a file of data acquisition integer data	*/
/*                                                                      */
/* To compile: cc -O -o deskip deskip.c  -lm				*/
/* To run: deskip 'sensitivity file' 'unscaled file' 'deskipped file'	*/
/* H.P. Gavin, Dept. of Civil Engineering, University of Michigan, 7-94 */
/************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../HPGnumlib/HPGutil.h"

#define HEAD_LINES	6	/* lines in the header of raw data file	*/

FILE    *fp_raw,		/* pointer to the raw data file		*/
	*fp_scl,		/* pointer to the re-sorted data file	*/
	*fp_sns;		/* pointer to the sensitivity data file	*/

int main ( int argc, char *argv[] )
{
	char    ch,             /* dummy character                      */
		line[MAXL],	/* a header line			*/
		title[MAXL];	/* title of the sensitivity data file	*/

	float   *r,		/* the raw data				*/
		sensi[16],	/* sensor sensitivities			*/
		sr,		/* sample rate ( samples / second )	*/
		range,		/* voltage range			*/
		*vector();     /* dynamic matrix memory allocation     */

	int     points=0,       /* number of data points                */
                nchan,		/* the number of channels               */
		start_chnl,	/* the first channel in the scan	*/
		stop_chnl,	/* the last channelin the scan		*/
		chn,		/* a channel number			*/
		len,		/* the length of a header line		*/
		skip_pt,	/* point where data acquisition skips	*/
		skip_ch,	/* number of channels that were skipped	*/
                i,j,k, q, sv;

	void	usr_inp(),	/* read input file, open output file	*/
        	free_vector();  /* dynamic matrix memory de-allocation  */


	usr_inp ( argc, argv, title, &sr, &range,
					&start_chnl, &stop_chnl, sensi );

                /*      Size  The  Data  And  The  Vectors      */

	nchan = stop_chnl - start_chnl + 1;

	for (i=1;i<=HEAD_LINES;i++)     /* scan through the header */
        	while (( ch = getc(fp_raw)) != '\n') ;
	while ( (ch=fscanf(fp_raw,"%f",&sr)) != EOF ) ++points;
	points /= nchan;
	printf (" %s: %d data points\n", argv[2], points);
	rewind (fp_raw);
	points *= nchan;

	r = vector( 1, points );

	for (i=1; i<=HEAD_LINES; i++) {		/* scan through the header */
		len = getLine ( fp_raw, MAXL , line);
		for (j=0; j<len; j++) putc ( line[j], fp_scl );
	}

        for ( k=1; k <= points; k++ )		/* read the data	*/
		ch=fscanf(fp_raw,"%f",&r[k]);

	fprintf(stderr," Last point before skip: ");
	sv=scanf("%d", &skip_pt);
	fprintf(stderr," Number of channels skipped: ");
	sv=scanf("%d", &skip_ch);

	k = 1;
	while ( k <= points ) {
		for (chn=start_chnl; chn <= stop_chnl; chn++) {
			fprintf(fp_scl,"%5.0f\t", r[k]);
			k++;
		}
		fprintf(fp_scl,"\n");
		if ( k == (skip_pt*nchan) + 1 ) {
			printf("skipping...\n");
			for (q=1; q <= skip_ch; q++)
				fprintf(fp_scl,"%5.0f\t", r[k+q-nchan-1]);
			for (chn=start_chnl+skip_ch; chn <= stop_chnl; chn++) {
				fprintf(fp_scl,"%5.0f\t", r[k]);
				k++;
			}
			fprintf(fp_scl,"\n");
		}
	}

	fclose (fp_sns);
	fclose (fp_raw);
	fclose (fp_scl);
	free_vector(r,1,points);
}

/*-----------------------------------------------------------------------------
USR_INP  -  read start up file, open input/output files			3feb94
		---- SENSITIVITY  DATA  FILE  FORMAT ----
  Descriptive Title (one line only)
  Sample rate (Hz):  Starting Channel:  Stopping Channel:
  Channel:		Sensitivity:		Units:
     "   :                   "     :              "  :
  Channel:		Sensitivity:		Units:
------------------------------------------------------------------------------*/

void usr_inp ( argc, argv, title, sr, range, start_chnl, stop_chnl, sensi )
int	argc;
char	*argv[];
float	*sr, *sensi, *range;
char	*title; 
int	*start_chnl, *stop_chnl;
{
	int	chn, nchan, i, sv;
	char	units[MAXL], line[MAXL];
	float	sns;

	if ((argc != 4) || ((fp_sns=fopen(argv[1],"r")) == NULL )) {
	    fprintf(stderr,"  usage: deskip [sens file] [infile] [outfile]\n ");
	    fprintf(stderr,"\t\t---- SENSITIVITY  DATA  FILE  FORMAT ----\n");
	    fprintf(stderr,"  Descriptive Title (one line only)\n");
	    fprintf(stderr,"  Sample rate (Hz):");
	    fprintf(stderr,"  Voltage range (V):");
	    fprintf(stderr,"  Starting Channel:");
	    fprintf(stderr,"  Stopping Channel:\n");
	    fprintf(stderr,"  Channel:\t\tSensitivity:\t\tUnits:\n");
	    fprintf(stderr,"    ''   :\t\t     ''    :\t\t  '' :\n");
	    fprintf(stderr,"  Channel:\t\tSensitivity:\t\tUnits:\n");
	    exit(1);
	}

	(void) getLine  ( fp_sns, MAXL , title);
	scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%f", sr );
	scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%f", range );
	scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%d", start_chnl );
	scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%d", stop_chnl );
	nchan = *stop_chnl - *start_chnl + 1;
	for (i=1; i <= nchan; i++) {
		scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%d", &chn );
		scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%f", &sns );
		scanLine ( fp_sns, MAXL, line, ':' );	sv=fscanf ( fp_sns, "%s", units );
		sensi[chn] = sns;
	}

        if (( fp_raw = fopen(argv[2], "r")) == NULL ) {
                printf("  error: cannot open file '%s'\n", argv[2]);
                printf("  usage: deskip <sensitivity file> <raw data file>");
                printf("  <deskipped file> \n");
                exit(2);
	}
        if (( fp_scl = fopen(argv[3], "w")) == NULL ) {
		printf("  error: cannot open file '%s'\n", argv[3]);
                printf("  usage: deskip <sensitivity file> <raw data file>\n");
                printf("  <deskipped file> \n");
                exit(3);
        }
	return;
}

