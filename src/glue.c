/* 
Program glue.c -  glues data file A to data file B, column-wise.

Data files A and B are multi-column data files, which have a set of 'header' 
lines starting with a # or % character.  Data file A can have more rows than
data file B, but can not have less rows than data file B.   

The 'glued' file will use the header of data file A and will place the columns
of data file B after the last column of data file A.  If data file A is longer
than data file B, the additional rows of data will be zero's.  

... currently, file B must be a one-column data file ... 

	to compile:	gcc -O -o glue glue.c 
	to run:		glue fileA fileB fileC

*/

#include <stdio.h>

#include "../../HPGnumlib/HPGutil.h"

int main ( int argc, char *argv[] )
{
	char	c='a', ch='a', 
		lineA[MAXL], 	/* a line from data file A */
		lineB[MAXL]; 	/* a line from data file B */

	FILE	*fpA,		/*  input file pointer - file A		*/
		*fpB,		/*  input file pointer - file B		*/
		*fpC;		/* output file pointer			*/

	int	i, j, len,
		head_lines_A=0, head_lines_B=0, rows_A=0, rows_B=0, rows;


	if ((fpA = fopen (argv[1], "r")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[1]);
		printf (" usage: glue file_A file_B outfile\n");
		exit(0);
	}
	if ((fpB = fopen (argv[2], "r")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[2]);
		printf (" usage: glue file_A file_B outfile\n");
		exit(0);
	}
	if ((fpC = fopen (argv[3], "w")) == NULL) {
		printf (" error: cannot open file '%s'\n", argv[3]);
		printf (" usage: glue file_A file_B outfile\n");
		exit(0);
	}

        head_lines_A = 0;      /* determine the size of the header in file A */
        do {
		ch = getc(fpA);
                if ( ch == '#' || ch == '%') {
                        ++head_lines_A;
                        while (( c = getc(fpA) ) != '\n' ) ;
                }
        } while ( ch == '#' || ch == '%' );
        rewind (fpA);

        head_lines_B = 0;      /* determine the size of the header  in file B */
        do {
		ch = getc(fpB);
                if ( ch == '#' || ch == '%' ) {
                        ++head_lines_B;
                        while (( c = getc(fpB) ) != '\n' ) ;
                }
        } while ( ch == '#' || ch == '%' );
        rewind (fpB);

	for (i=1;i<=head_lines_A;i++)		/* scan through file A */
                while (( ch = getc(fpA)) != '\n') ;

	for (i=1;i<=head_lines_B;i++)		/* scan through file B */
                while (( ch = getc(fpB)) != '\n') ;

	while ( getLine (fpA, MAXL, lineA) > 0 )	rows_A++;
	while ( getLine (fpB, MAXL, lineB) > 0 )	rows_B++;

	rewind(fpA);
	rewind(fpB);

	fprintf(fpC,"%% ==================================================\n");
	fprintf(fpC,"%% file %s\n", argv[1] );
        for (i=1; i<=head_lines_A; i++) { /* copy header A to file C */
                len = getLine ( fpA, MAXL , lineA);
		fprintf(fpC,"%s", lineA);
        }

	fprintf(fpC,"%% ==================================================\n");
	fprintf(fpC,"%% file %s\n", argv[2] );
	for (i=1;i<=head_lines_B;i++) {	 /* copy header B to file C */
                len = getLine ( fpB, MAXL, lineB );
		fprintf(fpC,"%s", lineB);
	}
	fprintf(fpC,"%% ==================================================\n");
	fprintf(fpC,"%% File '%s'   glued to file  '%s'\n", argv[2], argv[1] );

	rows = (rows_A < rows_B) ? rows_A : rows_B;

	for (i=1; i<=rows; i++) {
		len = getLine (fpA, MAXL, lineA);   
		len = getLine (fpB, MAXL, lineB);   
		fprintf(fpC,"%s\t%s", lineA, lineB);
	}

	/*
	for (i=rows_B+1; i<=rows_A; i++) {
		len = getLine (fpA, MAXL, lineA);   
		fprintf(fpC,"%s\t0\n", lineA );
	}
	*/

	printf(" rows in file A = %d   rows in file B = %d \n", rows_A, rows_B);

	fclose (fpA);
	fclose (fpB);
	fclose (fpC);
}



