/*******************************************************************************
Program xfer.c - calculates frequency response functions H1, H2, and Hv, and the
associated power spectra Gu, Gy, and Guy using data sets x & y. 
Averaging and Windowing are Incorporated.	 
  ==> Press, William H. et al.,  Numerical Recipes in C, 1st edition,
	(Cambridge, England: Cambridge University Press, 1988): Chapter 12.
  ==> Ramsey, K.A., ``Effective Measurements for Structural Dynamics Testing'',
	Sound and Vibration Magazine, (Nov. 1975): 24--35.
  ==> Rocklin, G.T., J. Crowley, and H. Vold, ``A comparison of H1, H2 and Hv
	Frequency Response Functions'', Proc. 3rd Int'l Modal Analysis Conf.,
	(Orlando, FL, Feb. 1985): 272--278.
  ==> Roth, P.R., ``How to Use the Spectrum and Coherence Function'',
	Sound and Vibration Magazine, (Jan. 1971): 10--14.
  ==> Vold, H., J. Crowley, and G.T. Rocklin, ``New Ways of Estimating Frequency
	Response Functions'', Sound and Vibration Magazine, (Nov. 1984): 34--38.

   to compile :	gcc -O -o xfer xfer.c fft.c complex.c NRutil.c -lm
   to run     :	xfer [data file]
	H.P. Gavin, Dept. of Civil Eng., Princeton University 		Feb 1991
*******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include "../../HPGnumlib/NRcomplex.h"
#include "../../HPGnumlib/NRutil.h"
#include "../../HPGnumlib/HPGsignal.h"
#include "../../HPGnumlib/NRcomplex.h"

float	head_lines;		/* number of lines in the file header	*/

FILE	*fp;		/* file pointer to input and output data files	*/

int main ( int argc, char *argv[] )
{
	char	ch;

	f_complex *Guy,		/* complex cross power spectrum		*/
		 *H;		/* comples frequency response function	*/

	float	*Gu, *Gy,	/* power spectra of x & y data files	*/
		*Coh,		/* coherence function			*/
		df;		/* frequency resolution interval	*/

	int	m = 512,	/* number of frequency values / segment	*/
		k = 5,		/* number of segments to be averaged	*/
		ovrlap = 1,	/* 1: overlap segments; 0: no overlap	*/
		points = 0,	/* number of data points ( 10,000 )	*/
		Ucol, Ycol,	/* columns for U and Y data	     */
		Ncol;		/* number of columns in the data file   */

	void	usr_inp(),      /* open input file, read control data   */
		spectra(),	/* returns spectra and cross spectra	*/
		frf(),		/* frequency response estimates:H1,H2,Hv*/
		save_data();	/* save the spectra to files		*/


	usr_inp ( argc, argv, &Ncol, &Ucol, &Ycol, &points, &df, &k, &m );

	Gu  =  vector(1,m);
	Gy  =  vector(1,m);
	Guy = Cvector(1,m);

	spectra ( Ncol, Ucol, Ycol, Gu,Gy,Guy, m,k, ovrlap );

	fclose (fp);

	H   = Cvector(1,m);
	Coh =  vector(1,m);

	frf ( Gu,Gy,Guy, H, Coh, df, m );

	save_data ( argc, argv, Ucol, Ycol, Gu, Gy, Guy, H, Coh, df, m );

	fprintf (stderr,"\n");

	free_vector(Gu,1,m);
	free_vector(Gy,1,m);
	free_Cvector(Guy,1,m);
	free_vector(Coh,1,m);
	free_Cvector(H,1,m);
}


/*------------------------------------------------------------------------------
USR_INP  -  open input file, read in control data			25jun93
------------------------------------------------------------------------------*/
void usr_inp ( argc, argv, Ncol, Ucol, Ycol, points, df, k, m )
int	argc;
char    *argv[];
float	*df;
int     *Ncol, *Ucol, *Ycol, *points, *k, *m;
{
	char    c = 'a', ch = 'a';
	float	data[20], sr;	/* vibration data and the sample rate	*/
	int	tg_pts,		/* target number of freq vals / segment	*/
		i, j = 1, sv;

	if ((fp = fopen(argv[1],"r")) == NULL ) {
	 fprintf(stderr,"  error: cannot open time series file '%s'\n", argv[1]);
	 fprintf(stderr,"  usage: xfer <uy file>");
	 fprintf(stderr," <N cols> <U-col> <Y-col> <sr> <df>");
	 fprintf(stderr," <G file> <H file>\n");
	 fprintf(stderr,"  ... where ... \n");
	 fprintf(stderr,"  <uy file> =  time series file name\n");
	 fprintf(stderr,"  <N cols>  =  number of colums in <uy file>\n");
	 fprintf(stderr,"  <U col>   =  input  data column in <uy file>\n");
	 fprintf(stderr,"  <Y col>   =  output data column in <uy file>\n");
	 fprintf(stderr,"  <sr>      =  sample rate of time series\n");
	 fprintf(stderr,"  <df>      =  desired frequency interval\n");
	 fprintf(stderr,"  <G file>  =  file name for spectra\n");
	 fprintf(stderr,"  <H file>  =  file name for frequency response function\n");
		exit(0);
	}

	if (argc < 3) {
		fprintf(stderr," Number of Columns: ");
		sv=scanf("%d", Ncol);
	} else	*Ncol = atoi(argv[2]);
	if (argc < 4) {
		fprintf(stderr," Input  data column (U): ");
		sv=scanf("%d", Ucol);
	} else	*Ucol = atoi(argv[3]);
	if (argc < 5) {
		fprintf(stderr," Output data column (Y): ");
		sv=scanf("%d", Ycol);
	} else	*Ycol = atoi(argv[4]);
	if (argc < 6) {
		fprintf (stderr," Sample rate (Hz): ");
		sv=scanf ("%f", &sr);
	} else	sr = atof(argv[5]);
	if (argc < 7) {
		fprintf (stderr," Target  frequency interval (Hz): "); 
		sv=scanf ("%f", df);
	} else	*df = atof(argv[6]);

	printf("  Ncol=%d  Ucol=%d  Ycol=%d sr=%.1f df=%.3f",
					*Ncol, *Ucol, *Ycol, sr, *df );

		/*	Size  The  Data		*/

	head_lines = 0;			/* determine the size of the header */
	do {
		ch = getc(fp);
		if ( ch == '#' || ch == '%' ) {
			++head_lines;
			while (( c = getc(fp) ) != '\n' ) ;
		}
	} while ( ch == '#' || ch == '%' );
	rewind (fp);

	for (i=1;i<=head_lines;i++)     /* scan through the header */
		while (( ch = getc(fp)) != '\n') ;

	do {
       		for (i=1; i<= *Ncol; i++) {
			j=fscanf(fp,"%f", &data[i] );
			if (j != 1) break;
		}
		if (j != 1) break;
		(*points)++;
	} while (j == 1);
	fprintf (stdout,"  %d data points\n", *points);
	rewind ( fp );

		/*	Choose  Frequency  Resolution	*/

	tg_pts = sr / *df;
	while ( *m < tg_pts) *m <<= 1; *m >>= 2;
	*k = *points / *m - 1;
	*df = sr / ( 2 * *m);

	fprintf (stdout,"  Computing %d averages of %d point spectra", *k, *m );
	fprintf (stdout," with %f Hz resolution ...\n", *df );
	if ( *k < 1 ) {
		fprintf (stdout,"%d is too small!\n", *k);
		exit(0);
	}
	return;
}


/*------------------------------------------------------------------------------
SPECTRA  -  Power Spectra Using The Fast Fourier Transform	
		Adapted from Numerical Recipes in C - Ch 12	
------------------------------------------------------------------------------*/
#define tpi	6.28318530717959

#define WINDOW(j,c,b) (0.5-0.5*cos(tpi*((j)-1)*(c)))	/* Hanning	*/
/* #define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))	/* Parzen	*/
/* #define WINDOW(j,a,b) (1.0-SQR((((j)-1)-(a))*(b)))	/* Welch	*/
/* #define WINDOW(j,a,b) 1.0	 			/* Square	*/

void spectra ( Ncol, Ucol, Ycol, Gu,Gy,Guy, m,k,ovrlap )
f_complex *Guy;
float	Gu[], Gy[];
int	m,k,ovrlap, Ncol, Ucol, Ycol;
{
	char	ch;
	int	n,tn,tn4,tn3, kk, i,j, tj, sv;	/* Useful Quantities	*/
	float	a,b,c, sumw=0.0, den=0.0,
		data[20],
		*x1, *x2, *y1, *y2,		/* Operating vectors	*/
		*Sx,*Sy;			/* FFT'd vectors	*/
	void segments ( int Ncol, int Ucol, int Ycol, float *x1, float * x2, float *y1, float *y2, float a, float b, float c, int  m, int ovrlap ); 

	n   = m+m;
	tn3 = (tn=n+n)+3;
	tn4 = tn3+1;
	a = m - 0.5;
	b = 1.0/(m+0.5);
	c = 1.0/(n-1.0);

	x1 = vector(1,m);
	x2 = vector(1,n);
	y1 = vector(1,m);
	y2 = vector(1,n);
	Sx = vector(1,tn);
	Sy = vector(1,tn);

	for (j=1; j<=n; j++) sumw += SQR(WINDOW(j,c,b));
	for (j=1; j<=m;  j++) Gu[j]=Gy[j]=Guy[j].r=Guy[j].i=0.0;

	for (i=1;i<=head_lines;i++)     /* scan through the header */
		while (( ch = getc(fp)) != '\n') ;

	if (ovrlap) {		/* Initialize the "save" half-buffer.	*/
		for (j=1; j<=m; j++) {
			for (i=1; i <= Ncol; i++) sv=fscanf(fp,"%f", &data[i] );
			x1[j] = data[Ucol];
			y1[j] = data[Ycol];
		}
	}

	for (kk=1; kk<=k; kk++) {

		segments ( Ncol, Ucol, Ycol, x1,x2,y1,y2, a,b,c, m, ovrlap );

				/* subtract a straight line from the data */
		deTrend ( x2, n );
		deTrend ( y2, n );

				/* Fourier transform the windowed data	*/
		twofft ( x2, y2, Sx, Sy, n );

			/* Sum the results into the previous segments	*/ 
		Gu[1] += (SQR(Sx[1]) + SQR(Sx[2]));
		Gy[1] += (SQR(Sy[1]) + SQR(Sy[2]));

		for (j=2; j<=m; j++) {
			tj = j+j;
			Gu[j] += (SQR(Sx[tj-1]) + SQR(Sx[tj])
				+ SQR(Sx[tn3-tj]) + SQR(Sx[tn4-tj]));
			Gy[j] += (SQR(Sy[tj-1]) + SQR(Sy[tj])
				+ SQR(Sy[tn3-tj]) + SQR(Sy[tn4-tj]));
			Guy[j].r += (Sx[tj-1]*Sy[tj-1] + Sx[tj]*Sy[tj]
				 	+ Sx[tn3-tj]*Sy[tn3-tj]
					+ Sx[tn4-tj]*Sy[tn4-tj]);
			Guy[j].i += (Sx[tj]*Sy[tj-1] - Sx[tj-1]*Sy[tj]
					- Sx[tn4-tj]*Sy[tn3-tj] 
					+ Sx[tn3-tj]*Sy[tn4-tj]);
		}
		den += sumw;
	}

	den *= tn;			/* Correct normalization.	*/
	for (j=1; j<=m; j++) {
		Gu[j] /= den;
		Gy[j] /= den;
		Guy[j].r /= den;
		Guy[j].i /= den;
	}
	free_vector(x1,1,m);
	free_vector(x2,1,n);
	free_vector(y1,1,m);
	free_vector(y2,1,n);
	free_vector(Sx,1,tn);
	free_vector(Sy,1,tn);
	return;
}


/*------------------------------------------------------------------------------
SEGMENTS  -  Gets two complete segments into the work space. 	See page 446
------------------------------------------------------------------------------*/
void segments ( int Ncol, int Ucol, int Ycol, float *x1, float * x2, float *y1, float *y2, float a, float b, float c, int  m, int ovrlap )
{
	float	data[20];
	int 	i,j, sv;

	if (ovrlap) {
		for (j=1; j<=m; j++) {
			x2[j]=x1[j];
			y2[j]=y1[j];
		}
		for (j=1; j<=m; j++) {
			for (i=1; i <= Ncol; i++)   sv=fscanf(fp,"%f", &data[i] );
			x1[j] = data[Ucol];
			y1[j] = data[Ycol];
		}
		for (j=1; j<=m; j++) {
			x2[m+j]=x1[j];
			y2[m+j]=y1[j];
		}
	} else {
		for (j=1; j<=m+m; j++) {
			for (i=1; i <= Ncol; i++)   sv=fscanf(fp,"%f", &data[i] );
			x2[j] = data[Ucol];
			y2[j] = data[Ycol];
		}
	}
	for (j=1; j<=m+m; j++) {	/* Apply the window to the data */
		x2[j] *= WINDOW(j,c,b);
		y2[j] *= WINDOW(j,c,b);
	}
	return;
}



/*------------------------------------------------------------------------------
FRF  -  Computes the frequency response function three ways using the
	auto power spectra Gu and Gy and the cross power spectrum Guy		
------------------------------------------------------------------------------*/
void frf ( Gu,Gy,Guy, H, Coh, df, m )
f_complex *Guy, *H;
float	*Gu,*Gy, *Coh, df;
int	m;
{
	char	base = 'n', frf_type = '1', disp = 'n';
	int	j, sv;
	f_complex H1, H2, Hv;
	float	 Guy2;

/*	printf(" Base - excited vibrations? (y/n): ");
	sv=scanf ( "%s", &base );
	printf(" Displacemnt frf from acceleration Data? (y/n): ");
	sv=scanf ( "%s", &disp );
	printf(" FRF estimate type: H1, Hv, or H2  (1,v,2): ");
	sv=scanf ( "%s", &frf_type );
*/
	
	for (j=1; j<=m; j++) {
						/* Coherence estimate	*/
		Guy2 = SQR ( Guy[j].r )  +  SQR ( Guy[j].i );
		Coh[j] = Guy2 / ( Gu[j] * Gy[j] );

						/* H1 estimate	*/
		H1.r =  Guy[j].r / Gu[j];
		H1.i = -Guy[j].i / Gu[j]; 
						/* H2 estimate	*/
		if (Guy2 == 0.0) H2.r = H2.i = 0.0;
		else {
			H2.r =  Gy[j]*Guy[j].r / Guy2;
			H2.i = -Gy[j]*Guy[j].i / Guy2;
		}
						/* Hv estimate	*/
		Hv   =  Csqrt ( Cmul ( H1 , H2 ) );
		if ( Guy[j].r * Hv.r < 0 )	Hv.r *= -1.0;
		if ( Guy[j].i * Hv.i > 0 )	Hv.i *= -1.0;

		if ( frf_type == '1' )		H[j] = H1;
		if ( frf_type == '2' )		H[j] = H2;
		if ( frf_type == 'v' )		H[j] = Hv;
	}

					/* subtract base acceleration */
	if (base == 'y')  for (j=1; j<=m; j++)	H[j].r -= 1.0;

	if (disp == 'y') {		/* Scale from Accel. to Disp.	*/
		for (j=2;j<=m;j++) {
			H[j].r /= ( 4.*PI*PI*(j-1)*(j-1)*df*df );
			H[j].i /= ( 4.*PI*PI*(j-1)*(j-1)*df*df );
		}
	}
	return;
}


/*-----------------------------------------------------------------------------
SAVE_DATA  -  save specra and frequency response function		10jan96
------------------------------------------------------------------------------*/
void save_data ( argc, argv, Ucol, Ycol, Gu, Gy, Guy, H, Coh, df, m )
int	argc;
char	*argv[];
f_complex *Guy, *H;
float	 *Gu, *Gy, *Coh, df;
int	 Ucol, Ycol, m;
{
	char	filename[64];
	int	j, sv;

	if (argc < 8) {
		fprintf (stderr," Power  Spectra   file  name: ");
		sv=scanf ("%s", filename); 
	} else strcpy(filename,argv[7]);


	if (( fp=fopen(filename,"w"))==NULL) {
		printf ("\n error: cannot open data file '%s'\n", filename);
		exit(0);
	}
	fprintf(fp,"%% power spectrum and cross spectrum file from data file");
	fprintf(fp," '%s'\n", argv[1]);
	fprintf(fp,"%% excitation column (U): %d,   response column (Y): %d,",
								Ucol, Ycol );
	fprintf(fp,"    delta-f= %f.3 Hz\n", df );
	fprintf(fp,"%% frequency    Gu	   Gy	   Guy-real");
	fprintf(fp,"     Guy-imag    coherence\n");
	for (j=1;j<=m;j++)
		fprintf (fp,"%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
			(j-1)*df,Gu[j],Gy[j],Guy[j].r,Guy[j].i,Coh[j]);
	fclose (fp);


	if (argc < 9) {
		fprintf (stderr," Transfer Function file name: ");
		sv=scanf ("%s", filename); 
	} else strcpy(filename,argv[8]);

	if (( fp=fopen(filename,"w"))==NULL) {
		printf ("\terror:cannot open data file '%s'\n",filename);
		exit(0);
	}
	fprintf(fp,"%% frequency response spectrum file from time history file");
	fprintf(fp," '%s'\n", argv[1]);
	fprintf(fp,"%% excitation column (U): %d,   response column (Y): %d,",
								Ucol, Ycol );
	fprintf(fp,"    delta-f= %f.3 Hz\n", df );
	fprintf(fp,"%% frequency\t  Hxy-real\t  Hxy-imag\n");
	for (j=2;j<=m;j++)
		fprintf (fp,"%11.4e %11.4e %11.4e\n", (j-1)*df,
/*					H[j].r, H[j].i);	/* real-imag */
			Cabs(H[j]), Carg(H[j]) * 180.0 / PI );	/* amp-phase */

	fclose (fp);

	return;
}
