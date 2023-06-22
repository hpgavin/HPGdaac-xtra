/******************************************************************************
 
Program chirpq.c - Generate a file of square waveform data with changing frequency and amplitude

  input:  sr      : sample rate (samples per second) 
          offset  : offset from zero
          ao      : initial amplitude 
          af      : final amplitude 
          fo      : initial frequency
          ff      : final frequncey
          T       : time duration of the data
          p       : amplitude distribution
          q       : frequency distribution
          outFile : output data file

  output: u       : random output signal

to compile: gcc -O -o chirpq chirpq.c HPGutil.c NRutil.c -lm 

to run:     chirpq sr offset ao af fo ff T p q outFile 

Henri Gavin, Dept. Civil Engineering, Duke University, 2022-12-04
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../HPGnumlib/HPGutil.h"
#include "../../HPGnumlib/NRutil.h"

#define PI       3.14159265358979323846264338327950288419716939937510
#define DAHI 0xFFFF
#define DALO 0x0000
#define DA00 0x0000
 
void display_parameters( FILE *fp, float sr, float offset, float ao, float af, float fo, float ff, float T,  float p, float q,  char *fn );

int main(int argc, char *argv[])
{
  char  fn[80],            // the output sweep filename
        usage[128];

  FILE *fp;                // file pointer to the output data file

  float  sr = 200.0,       // the sampe rate    
         offset = 0,     // static offset "set point"  
         ao = 50,          // starting amplitude
         af = 50,          // stopping amplitude  
         fo =  1.0,        // starting frequency limit
         ff = 10.0,        // stopping frequency limit
         T = 30.0,         // the time to sweep the frequency band
         p = 1.0,          // power of geometric increase in freq.
         q = 1.0, r = 1.0, // constants in amplitude envelope
         time  = 0.0,      // the current global time
         cycles,           // total number of cycles
         ai = 0.0,         // the current amplitude
         phase,            
        *u;                // vector of sweep data

  int    nTu, nTd,         // points in the taper stages
         nTs,              // number of points in sweep
         nT,               // total number of points in the record
         t0 = 0,           // the time at the start of the chirp
         t = 0,            // index of the data record
         sign = 0;         // the instantaneous sign 

  void   exit();

  sprintf(usage," chirpq sr offset ao af fo ff T p q outFile");

  printf("\n  to run:  %s \n\n", usage );

  // default values 
  if (argc < 11) sprintf(fn,"chirpq.dat"); else strcpy( fn, argv[10]);
  if (argc < 10) q = 1;                   else q      = atof( argv[9] );
  if (argc <  9) p = 1;                   else p      = atof( argv[8] );
  if (argc <  8) T  = 30.0;               else T      = atof( argv[7] );
  if (argc <  7) ff = 10.0;               else ff     = atof( argv[6] );
  if (argc <  6) fo =  1.0;               else fo     = atof( argv[5] );
  if (argc <  5) af = 20.0;               else af     = atof( argv[4] );
  if (argc <  4) ao = 50.0;               else ao     = atof( argv[3] );
  if (argc <  3) offset = 50.0;           else offset = atof( argv[2] );
  if (argc <  2) sr = 200.0;              else sr     = atof( argv[1] );

  if (ff/sr > 0.051 || fo/sr > 0.051) // want more than 20 points / cycle 
    printf("  warning: sample rate is too low for hydraulic actuators!\n");

  fp = openFile("./",fn,"w",usage);

  display_parameters( stdout , sr, offset, ao, af, fo, ff, T,  p, q, fn );
  display_parameters( fp     , sr, offset, ao, af, fo, ff, T,  p, q, fn );
 
  ao     = ao/100.0;              // initial amplitude - fraction of full scale
  af     = af/100.0;              //  final  amplitude - fraction of full scale
  offset = offset/100.0;          //      offset       - fraction of full scale

  nT  = (int)(T*sr);              // total number of points 
  t0  = (int)(1.00*sr);           // number of points for offset adjustment
  nTu = (int)(0.10*nT);           // number of points to taper up
  nTd = (int)(0.10*nT);           // number of points to taper down
  nTs = (int)(nT-2*t0);           // number of points in the sine sweep
  cycles = nTs * ( fo + (ff-fo)/(p+1));  // number of cycles in the record

  p = fabs(p);
  if (fo<ff && p < 1) p = 1/p;
  if (fo>ff && p > 1) p = 1/p;

  r   =  log(ao/af) * pow(T, -q);   // amplitude exponential in time 

  u = vector ( 1 , nT );            // allocate memory to vector u

  for ( t = 1; t <= nT; t++ ) u[t] = DA00;

  // frequency sweep via instantaneous phase angle
 
  sign = +1;
  for ( t=t0 ; t <= t0 + nTs; t++ ) {
    ai = ao + (af - ao) * pow(((float)(t)/(float)(cycles)),q);
    phase = 2*PI * ((t-t0)*fo/sr + ( ff - fo )*pow((t-t0)/sr,p+1) /
                                   ( (p+1)*pow(T,p) ));

    if ( cos(phase) > 0 ) sign = +1;
    if ( cos(phase) < 0 ) sign = -1;

    u[t] = sign*ai;
  }
  
  // taper up the inital section and taper down the final section
  for ( t=0; t<=nTu; t++ ) u[t0+t]    *= 0.5*( 1.0-cos(PI*t/(nTu+1)) );
  for ( t=0; t<=nTd; t++ ) u[nT-t0-t] *= 0.5*( 1.0-cos(PI*t/(nTd+1)) );

  // scale to least significant bits, LSB
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] + offset) * (DAHI - DALO);

  // offset adjustment
  for ( t=0; t<=t0; t++ )  u[ 1+t] *= 0.5*( 1.0-cos(PI*t/(t0+1)) );
  for ( t=0; t<=t0; t++ )  u[nT-t] *= 0.5*( 1.0-cos(PI*t/(t0+1)) );
  
  // clip to the range DALO to DAHI 
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] > DAHI ) ? DAHI : u[t];
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] < DALO ) ? DALO : u[t];

  // write to the output file
  for ( t=1; t<=nT; t++ )  fprintf( fp,"%7d\n", (int)(round(u[t])) );

  free_vector( u , 1 , nT );

  fclose (fp);

  return(0);
}

void display_parameters( FILE *fp, float sr, float offset, float ao, float af, float fo, float ff, float T,  float p, float q,  char *fn )
{
  fprintf(fp,"%%   sr      : sampe rate \t\t(samples/sec)\t[ %4.0f ]\n", sr );
  fprintf(fp,"%%   offset  : offset \t\t\t\t(%%FS)\t[ %4.1f ]\n", offset );
  fprintf(fp,"%%   ao      : starting ampliutde \t\t(%%FS)\t[ %4.1f ]\n", ao );
  fprintf(fp,"%%   af      : stopping amplitude \t\t(%%FS)\t[ %4.1f ]\n", af );
  fprintf(fp,"%%   fo      : starting frequency \t\t(Hz)\t[ %4.1f ]\n", fo );
  fprintf(fp,"%%   ff      : stopping frequency \t\t(Hz)\t[ %4.1f ]\n", ff );
  fprintf(fp,"%%   T       : time duration \t\t\t(sec)\t[ %4.1f ]\n", T );
  fprintf(fp,"%%   p       : power of geometric frequency sweep \t[ %4.1f ]\n", p );
  fprintf(fp,"%%   q       : time constant for amplitude variation \t[ %4.1f ]\n", q );
  fprintf(fp,"%%   outFile : name of the output file\t\t\t[ %s ]\n\n", fn);

  return;
}
   
