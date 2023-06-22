/*****************************************************************************

Program ktrand.c - synthetic earthquake ground motion with the Kannai Tajimi spectrum

  input:  sr      : sample rate,         samples per second
          offset  : static offset                       %FS
          RMS     : rms amplitude ( > 0 )               %FS
          fg      : natural frequency of ground motion,  Hz
          zg      : damping ratio of ground motion
          aa      : envelope rise  parameter
          tt      : envelope decay parameter
          T       : time duration,                        s
          seed    : seed of the random number generator
          file    : name of output file

  output: Accel   : ground acceleration

  to compile:  gcc -O -o ktrand ktrand.c NRutil.c HPGrandom.c  -lm
      to run:  ktrand  sr RMS fg zg aa tt T seed datafile
 
  Henri Gavin, Civil and Environ. Engineering, Duke University, March 2007
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../HPGnumlib/NRutil.h"
#include "../../HPGnumlib/HPGutil.h"
#include "../../HPGnumlib/HPGsignal.h"
#include "../../HPGnumlib/HPGrandom.h"
#include "../../HPGnumlib/HPGcholesky.h"

#define NF 512	/* number of discrete frequencies in spectrum	*/
#define PI 3.14159265358979323846264338327950288419716939937510

#define DAHI 0xFFFF
#define DALO 0x0000
#define DA00 0x0000


void display_parameters( FILE *fp, float sr, float offset, float RMS, float fg, float zg, float aa, float tt, float T, int seed, char *fn );


int main(int argc, char *argv[])
{

  char	fn[80],		// the output random data filename
	  ch;

  FILE	*fp;		// file pointer to the output data file	

  float	sr, dt,	        // sample rate, time step, s
        a, v, d, Dd,    // accel, veloc, displ
  	Ms,Cs,Ks,Bs,    // structural mass, damping, stiffness, input
        Mo,Co,Ko,       // linear acceleration coefficients 	
	fg, zg, 	// frequency and damping of ground
        aa, tt,         // envelope parameters
	RMS,		// root mean square of pre-envelope process
        T,              // duration
        offset,         // offset 
        ut0, ut1,       // unit Gaussean noise
        mAV,            // max absolute value of terms in a vector 
        *u;             // vector of output process

  int	seed,		// seed for the random number generator
        idx, 
	t, nT, t0, nTu, nTd; 

  char	usage[80];


  sprintf(usage," ktrand sr offset fg zg aa tt T seed file_name");

  printf("\n  to run:  %s \n\n", usage );

//	 Default Values 
  if (argc < 11) strcpy(fn,"ktrand.dat"); else strcpy( fn, argv[10] );
  if (argc < 10) seed = time(0);	  else seed    = atoi( argv[9] );
  if (argc < 9) T = 30.0;           	  else T       = atof( argv[8] );
  if (argc < 8) tt = 4.0;          	  else tt      = atof( argv[7] );
  if (argc < 7) aa = 2.0;          	  else aa      = atof( argv[6] );
  if (argc < 6) zg = 0.9;          	  else zg      = atof( argv[5] );
  if (argc < 5) fg =  1.5;	          else fg      = atof( argv[4] );
  if (argc < 4) RMS = 0.23;	          else RMS     = atof( argv[3] );
  if (argc < 3) offset = 50.0;	          else offset  = atof( argv[2] );
  if (argc < 2) sr = 100.0;	          else sr      = atof( argv[1] );

  fp = openFile("./",fn,"w",usage);

  display_parameters( stdout , sr, offset, RMS, fg, zg, aa, tt, T, seed, fn );
  display_parameters( fp     , sr, offset, RMS, fg, zg, aa, tt, T, seed, fn );

  t0  = (int)(1.0*sr);                // number of points for offset adjustment
  nTu = (int)(0.1*nT);                // number of points to taper up
  nTd = (int)(0.1*nT);                // number of points to taper down

  srand48(seed);

  dt = 1.0 / sr; 

  nT = floor( sr * T );

  u  = vector(1,nT);

  Ms = 1.0;
  Cs = 4.0*PI*zg*fg;
  Ks = 4.0*PI*PI*fg*fg;
  Bs = RMS / sqrt(2.0*PI*zg*fg); 

  //  linear acceleration method
  Mo = 3.0*Ms       +     Cs*dt/2.0;
  Co = 6.0*Ms/dt    + 3.0*Cs;
  Ko = 6.0*Ms/dt/dt + 3.0*Cs/dt     + Ks;

  //  initial conditions 
  a = 0.0;
  v = 0.0;
  d = 0.0;

  for (t=1; t<=nT; t++) u[t] = DA00;  //  initialize output vector

  ut0 = 0.0;
  for ( t=t0; t<nT-t0; t++ ) {
    ut1  = ut0;
    ut0  = invnorm( drand48() ) / sqrt(dt);
    Dd   = ( Mo*a + Co*v + Bs*(ut0-ut1) ) / Ko;
    d   += Dd;
    v    = (-2.0*v - a*dt/2.0 + 3.0*Dd/dt);
    a    = ( Bs*ut0 - Cs*v - Ks*d ) / Ms ;    // e.o.m.
    u[t] = 4.0*PI*fg*zg*v; 
  }

  // envelope
  for ( t=t0; t<nT-t0; t++ ) 
    u[t] *= pow( (t-t0)*dt/(aa*tt) , aa ) * exp( aa - (t-t0)*dt/tt );

  // taper up the inital section and taper down the final section
  for ( t=0; t<=nTu; t++ ) u[t0+t]    *= 0.5*( 1.0-cos(PI*t/(nTu+1)) );
  for ( t=0; t<=nTd; t++ ) u[nT-t0-t] *= 0.5*( 1.0-cos(PI*t/(nTd+1)) );
  
  // scale to least significant bits, LSB
  mAV = maxAbsV( u, nT, &idx );
  offset /= 100;  // convert from percentage to fraction
  for ( t=1; t<=nT; t++ )  u[t] = (u[t]*0.5/mAV + offset) * (DAHI - DALO);

  // smooth offset adustment
  for ( t=0; t<=t0; t++ )  u[ 1+t] *= 0.5*( 1.0-cos(PI*t/(t0+1)) );
  for ( t=0; t<=t0; t++ )  u[nT-t] *= 0.5*( 1.0-cos(PI*t/(t0+1)) );

  // clip to the range DALO to DAHI
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] > DAHI ) ? DAHI : u[t];
  for ( t=1; t<=nT; t++ )  u[t] = (u[t] < DALO ) ? DALO : u[t];
  
  // write to the output file
  for ( t=1; t<=nT; t++ )  fprintf(fp,"%7d\n", (int)(round(u[t])) );

  fclose (fp);

  free_vector(u,1,nT);

  return(0);
}


void display_parameters( FILE *fp, float sr, float offset, float RMS, float fg, float zg, float aa, float tt, float T, int seed, char *fn )
{
  fprintf(fp,"%%     sr      : sample rate \t\t\t\t[%.4f]\n", sr );
  fprintf(fp,"%%  offset     : offset shift \t\t\t\t[%.4f]\n", offset );
  fprintf(fp,"%%     RMS     : root mean square of motion\t\t[%3.1f]\n", RMS );
  fprintf(fp,"%%     fg      : natural frequency of ground motion\t[%3.1f]\n", fg );
  fprintf(fp,"%%     zg      : damping ratio of ground motion\t\t[%3.1f]\n", zg);
  fprintf(fp,"%%     aa      : envelope rise  time parameter \t\t[%3.1f]\n", aa);
  fprintf(fp,"%%     tt      : envelope decay time parameter \t\t[%3.1f]\n", tt);
  fprintf(fp,"%%     T       : time duration \t\t\t\t[%f]\n",T);
  fprintf(fp,"%%     seed    : seed of the random number generator\t[%d]\n",seed);
  fprintf(fp,"%%     outfile : name of the output file\t\t\t[ktrand.out]\n\n");
}

