/*****************************************************************************

Program bwrand.c - Generate a vector of gaussean band-limited random noise 

  input:  sr      : sample rate (samples per second) 
          offset  : offset from zero
          ampl    : ampliutde (percent full scale)
          f_lo    : low  frequency limit of band-limited noise
          f_hi    : high frequency limit of band-limited noise
          T       : time duration of the data
          Seed    : seed of the random number generator
          outFile : output data file

  output: u       : random output signal

to compile: 

gcc -O -o bwrand bwrand.c NRutil.c HPGutil.c HPGsignal.c HPGmatrix.c -lm

to run  :

bwrand  sr offset  f_lo  f_hi  T   seed   outfile
 
Henri Gavin, Civil and Environ. Engineering, Duke Univsersity, March 2007
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../HPGnumlib/HPGutil.h"
#include "../../HPGnumlib/HPGsignal.h"
#include "../../HPGnumlib/HPGmatrix.h"
#include "../../HPGnumlib/NRutil.h"

#define DAHI 0xFFFF
#define DALO 0x0000
#define DA00 0x0000

void display_parameters( FILE *fp, float sr, float offset,float ampl, float f_lo,float f_hi,float fw, float T, int seed, char *fn );

int main ( int argc,char *argv[] )
{

  char  fn[80],             // the output random data filename  
        usage[128];         //  how to run the program 

  FILE  *fp;                // file pointer to the output data file  
 
  float  sr = 200.0,        // sample rate, samples per second  
         offset = 0.0,      //  offset                              
         ampl = 50.0,       // amplitude (%FS)                      
         f_lo, f_hi,        // frequency band width      
         fw = 0,            // frequency weighting +1, 0 , -1
         T,                 // duration of the record     
         df,                // frequency increment      
         phase,             // random phase anlges      
         freq,              // frequency values      
         *u,                // output time series
         mAV = 0;
//       min,max,avg,rms,t_min,t_max; // statistics of a record  

  int    seed,              // seed for the random number generator  
         nTu, nTd,          // pointes to taper up and taper down 
         t, t0, nT,         // time indices
         q, q_lo, q_hi, nF, // frequency indices
         idx;               // index of a maximum absolute value  

  sprintf(usage," bwrand sr offset ampl f_lo f_hi fw T seed file_name");

  printf("\n  to run:  %s \n\n", usage );

  // default parameters
  if (argc <10) sprintf(fn,"bwrand.dat"); else strcpy( fn, argv[9]);
  if (argc < 9) seed = time(0);           else seed   = atoi( argv[8] );
  if (argc < 8) T = 60;                   else T      = atof( argv[7] );
  if (argc < 7) fw = 0;                   else fw     = atof( argv[6] );
  if (argc < 6) f_hi = 020.0;             else f_hi   = atof( argv[5] );
  if (argc < 5) f_lo = 0.5;               else f_lo   = atof( argv[4] );
  if (argc < 4) ampl = 50;                else ampl   = atof( argv[3] );
  if (argc < 3) offset = 0.0;             else offset = atof( argv[2] );
  if (argc < 2) sr = 200.0;               else sr     = atof( argv[1] );
 
  if ( f_hi/sr > 0.051)  // we need more than 20 points / cycle 
    printf("  warning: sample rate is too low for hydraulic actuators!\n");

  fp = openFile("./",fn,"w",usage);
  
  display_parameters( stdout , sr, offset,ampl, f_lo,f_hi,fw, T, seed, fn );
  display_parameters( fp     , sr, offset,ampl, f_lo,f_hi,fw, T, seed, fn );

  srand48(seed);                      // seed the random number generator

  ampl   = ampl/100.0;                // amplitude - fraction of full scale
  offset = offset/100.0;              //   offset  - fraction of full scale

  nT = (int)(T*sr);                   // total number of points 
  nF = floor(nT/2.0);                 // number of frequency indices
  df = sr/nT;                         // frequency increment

  t0  = (int)(1.0*sr);                // number of points for offset adjustment
  nTu = (int)(0.1*nT);                // number of points to taper up
  nTd = (int)(0.1*nT);                // number of points to taper down

  q_lo = (int)(round(f_lo/df));       // index of low  frequency
  q_hi = (int)(round(f_hi/df));       // index of high frequency

  u = vector( 1 , nT );               // allocate memory to vector u

  for (t=1; t<=nT; t++) u[t] = DA00;  // initialize output vector

  // random noise via Fourier series with random phases
  for ( q=q_lo; q<=q_hi; q++ ) {
    freq  = q*df;                // frequency, Hz
    phase = 2*PI*drand48();      // uniformly distributed random phase, rad
    for ( t=t0; t<=nT-t0; t++ ) {
      u[t] += sqrt(2.0/nF) * pow(2*PI*freq,fw) * cos( 2*PI*freq*t/sr + phase );
    }
  }

  // taper up the inital section and taper down the final section
  for ( t=0; t<=nTu; t++ ) u[t0+t]    *= 0.5*( 1.0-cos(PI*t/(nTu+1)) );
  for ( t=0; t<=nTd; t++ ) u[nT-t0-t] *= 0.5*( 1.0-cos(PI*t/(nTd+1)) );

  // scale to least significant bits, LSB
  mAV = maxAbsV( u, nT, &idx );
  for ( t=1; t<=nT; t++ )  u[t] = (u[t]*ampl/mAV + offset) * (DAHI - DALO);

  // smooth offset adustment
  for ( t=0; t<=t0; t++ )  u[ 1+t] *= 0.5*( 1.0-cos(PI*t/(t0+1)) );
  for ( t=0; t<=t0; t++ )  u[nT-t] *= 0.5*( 1.0-cos(PI*t/(t0+1)) );

  // clip to the range DALO to DAHI
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] > DAHI ) ? DAHI : u[t];
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] < DALO ) ? DALO : u[t];

  // write to the output file
  for ( t=1; t<=nT; t++ )  fprintf(fp,"%7d\n", (int)(round(u[t])) );

  fclose (fp);

  free_vector ( u , 1 , nT );

  return(0);
}


void display_parameters( FILE *fp, float sr, float offset,float ampl, float f_lo,float f_hi,float fw, float T, int seed, char *fn )
{
  fprintf(fp,"%%   sr      : sample rate\t\t(samples/sec)\t[ %5.1f ]\n", sr );
  fprintf(fp,"%%   offset  : offset\t\t\t\t(%%FS)\t[ %5.1f ]\n", offset );
  fprintf(fp,"%%   ampl    : amplitude\t\t\t\t(%%fs)\t[ %5.1f ]\n", ampl );
  fprintf(fp,"%%   f_lo    : low  frequency in band\t\t(Hz)\t[ %5.1f ]\n", f_lo );
  fprintf(fp,"%%   f_hi    : high frequency in band\t\t(Hz)\t[ %5.1f ]\n", f_hi );
  fprintf(fp,"%%   fw      : frequency weighting \t-2 < fw < +2\t[ %5.1f ]\n", fw );
  fprintf(fp,"%%   T       : duration of the record\t\t(sec)\t[ %5.1f ]\n", T );
  fprintf(fp,"%%   seed    : seed of the random number generator\t[ %5d ]\n",seed);
  fprintf(fp,"%%   outfile : name of the output file\t\t\t[ %s ]\n\n",fn);

  return;
}
