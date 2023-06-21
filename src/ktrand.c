/*****************************************************************************

Program ktrand.c - synthetic earthquake ground motion with the Kannai Tajimi spectrum

  input:  Sv      : peak amplitude ( > 0 )
          Fg      : natural frequency of ground motion
          Zg      : damping ratio of ground motion
          Vp      : pulse amplitude ( > 0 )
          Tp      : pulse period
          Nc      : pulse cycles
          points  : number of points in the vector  ( > 0 )
          delta_t : time step interval  
          seed    : seed of the random number generator
          file    : name of output file

  output: Accel   : ground acceleration
          Veloc   : ground velocity
          Displ   : ground displacement

  to compile:  gcc -O -o ktrand ktrand.c NRutil.c HPGsignal.c fft.c -lm
      to run:  ktrand  Sv Fg Zg Vp Tp Nc points delta_t seed file
 
  Henri Gavin, Civil and Environ. Engineering, Duke University, March 2007
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../HPGnumlib/NRutil.h"
#include "../../HPGnumlib/HPGsignal.h"
#include "../../HPGnumlib/HPGmatrix.h"

#define NF 512	/* number of discrete frequencies in spectrum	*/
#define PI 3.14159265358979323846264338327950288419716939937510


int main(int argc, char *argv[])
{

char	fn[80],		/* the output random data filename	*/
	ch;

FILE	*fp;		/* file pointer to the output data file	*/

float	phase[NF],	/* random phase anlges			*/
	freq[NF],	/* frequencies				*/
	ampl[NF],	/* amplitudes				*/
	delta_f,	/* frequency increment			*/
	sr, delta_t,	/* sample rate				*/
	Sv = 50.0,	/* spectral velocity			*/
	fg=2.5, zg=0.6,	/* frequency and damping of ground	*/
	Vp, Tp, Nc, ts,	/* pulse amplitude,period,cycles,time	*/
	f_min, f_max,	/* frequency bandwidth of simulation	*/
	P1, P2,	Pc,	/* points for ramping up and ramping down */
        T,              // duration
        offset,         // offset 
	*accel, *veloc, *displ,	/* acceleration, velocity and displacement */
	scaleA,		/* peak accel amplitude scaling		*/
	scaleV,		/* peak veloc amplitude scaling		*/
	scaleD,		/* peak displ amplitude scaling		*/
	scaleT = 1.0,	/* time scale factor for shake table	*/
	min,max,avg,rms,t_min,t_max; /* statistics of a record	*/

int	seed,		/* seed for the random number generator	*/ 
	idx,		/* index of a maximum absolute value	*/
	f, t, points;	/* number of points in the output file	*/

char	usage[80];


  sprintf(usage," ktrand sr offset fg zg T seed file_name");

//	 Default Values 
  if (argc < 9) strcpy(fn,"ktrand.dat"); else strcpy( fn, argv[8] );
  if (argc < 8) seed = time(0);	        else seed    = atoi( argv[7] );
  if (argc < 7) T = 30.0;         	else T       = atof( argv[6] );
  if (argc < 5) zg = 0.6;        	else zg      = atof( argv[5] );
  if (argc < 4) fg = 100;	        else Vp	     = atof( argv[4] );
  if (argc < 3) offset = 00.0;	        else offset  = atof( argv[2] );
  if (argc < 2) sr = 50;	        else sr      = atof( argv[1] );

  printf("  to run:  %s\n\n", usage );
  printf("     sr      : sample rate \t\t\t[%.4f]\n", sr );
  printf("  offset     : sample rate \t\t\t[%.4f]\n", sr );
  printf("     Fg      : natural frequency of ground motion\t[%3.1f]\n", fg );
  printf("     Zg      : damping ratio of ground motion\t\t[%3.1f]\n", zg);
  printf("     T       : time duration in the vector\t\t[%d]\n",points);
  printf("     seed    : seed of the random number generator\t[%d]\n",seed);
  printf("     outfile : name of the output file\t\t\t[ktrand.out]\n\n");

  delta_t = 1.0 / sr; 

  srand48(seed);
  accel = vector(1,points);
  veloc = vector(1,points);
  displ = vector(1,points);

  if (( fp = fopen(fn, "w" )) == NULL ) {
    printf ("  error: cannot open output file '%s' \n", fn );
    exit(1);
  }
	
	P1 =  ceil(0.10*points);	
	P2 = floor(0.50*points);
	Pc = 0.10*points;

	f_min = scaleT * 10.0 * sr / points;
	f_min = ( f_min > fg/10.0 ) ? fg/10.0 : f_min;
	f_max = 0.2 * sr;

	delta_f = (f_max - f_min)/NF;

	for (f=0;f<NF;f++) {
		phase[f] = 2*PI*drand48();	
		freq[f]  = f_min + f*delta_f;
 		ampl[f]  = sqrt( ( 1.0 + pow(2*zg*freq[f]/fg,2.) ) / 
		( pow(1-pow(freq[f]/fg,2.),2.) + pow(2*zg*freq[f]/fg,2.) ) );

	}

	for (t=1;t<=points;t++) {
	    accel[t] = 0;
	    for (f=0;f<NF;f++) 
		  accel[t] += ampl[f]*sin(2*PI*freq[f]*t*delta_t + phase[f]);
	    accel[t] *= sqrt(2.0/NF);
	    if (t<P1) accel[t] *= pow( ((t-1)/P1), 2.0);
	    if (t>P2) accel[t] *= exp(-(t-P2)/Pc);
	}

	veloc = baseLine( cumTick( accel, points, delta_t ),points,1,points,10);
	scaleV = Sv / maxAbsV(veloc,points, &idx);

	max = 0.0;
	for (t=1;t<=points;t++) {		/* add pulse	*/
		veloc[t] *= scaleV;
		ts = (t - 0.2*points)*delta_t*2*PI/Tp/Nc;  /* scaled time */
		veloc[t] += Vp*exp(-ts*ts/15) * (cos(Nc*ts)); /* even pulse */
//		veloc[t] += Vp*exp(-ts*ts/15) * (sin(Nc*ts)); /* even pulse */
		if ( fabs(veloc[t]) > max ) {
			max = abs(veloc[t]);
			if ( veloc[t] > 0 )	rms =  1.0;
			else			rms = -1.0;
		}
	}
	for (t=1;t<=points;t++) veloc[t] *= rms;

	accel = cDiff( veloc, points, delta_t );
	displ = baseLine( cumTick( veloc, points, delta_t ),points,1,points,10);
	
					/* print Header Information */
	fprintf(fp,"%%   Simulated Earthquake Time History : ");
	fprintf(fp,"    Kannai-Tajimi Spectrum\n");
	fprintf(fp,"%%   High-Freq Veloc:   %8.4f \n",Sv);
	fprintf(fp,"%%   Ground frequency:  %8.4f Hertz (prototype scale)\n",fg);
	fprintf(fp,"%%   Ground damping:    %8.4f percent \n", zg*100 );
	fprintf(fp,"%%   Pulse Velocity:    %8.4f \n", Vp );
	fprintf(fp,"%%   Pulse Period:      %8.4f \n", Tp );
	fprintf(fp,"%%   Pulse Cycles:      %8.4f \n", Nc );
	fprintf(fp,"%%   Number of points:  %d \n", points );
	fprintf(fp,"%%   Time Step:         %8.4f seconds \n", delta_t );
	fprintf(fp,"%%   seed (int):       %d \n", seed );
	fprintf(fp,"%%   Time Scale Factor: %8.4f \n", scaleT );
	fprintf(fp,"%%   Minimum frequency: %8.4f Hertz (model scale) \n",f_min);
	fprintf(fp,"%%   Maximum frequency: %8.4f Hertz (model scale) \n",f_max);
	fprintf(fp,"%%   Output filename:  %s \n", fn );


/*
	fprintf(fp,"%% _________________________________________________________\n");
	fprintf(fp,"%%         Minimum  Maximum  Average    RMS    Tmin    Tmax\n");
	limits(accel, points, &min, &max, &avg, &rms, &t_min, &t_max, delta_t );
	fprintf(fp,"%% accel: %8.3f %8.3f  %7.4f %7.3f %7.4f %7.4f\n",
					min, max, avg, rms, t_min, t_max );
	limits(veloc, points, &min, &max, &avg, &rms, &t_min, &t_max, delta_t );
	fprintf(fp,"%% veloc: %8.3f %8.3f  %7.4f %7.3f %7.4f %7.4f\n",
					min, max, avg, rms, t_min, t_max );
	limits(displ, points, &min, &max, &avg, &rms, &t_min, &t_max, delta_t );
	fprintf(fp,"%% displ: %8.3f %8.3f  %7.4f %7.3f %7.4f %7.4f\n",
					min, max, avg, rms, t_min, t_max );
	fprintf(fp,"%% _________________________________________________________\n");
	fprintf(fp,"%%\n");
*/

					/* print accel, veloc, displ */
	fprintf(fp,"%%  ACCEL cm/s/s    VELOC cm/s      DISPL cm \n");
	fprintf(fp,"%%\n");
 	for (t=1; t<=points; t++)		
		fprintf(fp,"%12.4e\t%12.4e\t%12.4e\n",
						accel[t], veloc[t], displ[t]);
 
					/* print displacement only */
/*
printf(" save acceleration (a)  velocity (v)  or   displacement (d)? " );
	scanf("%s", &ch);
	if (ch == 'a') {
		fprintf(fp,"%% acceleration time record ... cm/sec/sec \n");
		scaleA = 2000.0 / maxAbsV(accel,points, &idx);
		for (t=1; t<=points; t++)
			fprintf(fp,"%8.0f\n", accel[t]*scaleA );
	} else if (ch == 'v') {
		fprintf(fp,"%% velocity time record ... cm/sec \n");
		for (t=1; t<=points; t++)
			fprintf(fp,"%10.0f\n", veloc[t]*scaleV );
	} else if (ch == 'd') {
		scaleD = 2000.0 / maxAbsV(displ,points, &idx);
		fprintf(fp,"%% displacement time record ...  \n");
		for (t=1; t<=points; t++)
			 fprintf(fp,"%8.0f\n", displ[t]*scaleD );
	}
*/

	printf(" max. accel. = %f \n", maxAbsV(accel,points, &idx));
	printf(" max. veloc. = %f \n", maxAbsV(veloc,points, &idx));
	printf(" max. displ. = %f \n", maxAbsV(displ,points, &idx));

	fclose (fp);

	free_vector(accel,1,points);
	free_vector(veloc,1,points);
	free_vector(displ,1,points);

	return(0);
}

