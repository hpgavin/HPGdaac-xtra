/******************************************************************************
pulse.c - Generate a file with pulse-like waveform data. 

 to compile: gcc -O -o pulse pulse.c HPGsignal.c HPGmatrix.c fft.c gamma.c NRutil.c -lm

  input:  sr      : sample rate, samples per second 
          Tp      : period of the pulse
          Nc      : number of cycles in the pulse
          T       : duration of the record
          Vp      : pulse amplitude ( > 0ber generator
          file    : name of output file

 to run:     pulse sr Tp Nc T pulse.dat 


Henri Gavin, Dept. Civil Engineering, Duke University, 18 Dec. 2012
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../../HPGnumlib/HPGsignal.h"
#include "../../HPGnumlib/NRutil.h"
#include "../../HPGnumlib/HPGmatrix.h"
// #include "gamma.h"

int main ( int argc, char *argv[] )
{
	char	pulse_fn[80];	/* the output pulse motion filename	*/

	FILE	*pulse_fp;	/* file pointer t0 the output data file	*/

	float	pulse_time = 0.0, /* the time t0 sweep the frequency band */
		t = 0.0,	/* the current global time		*/
		t0 = 0.50,	/* the time at the start of the pulse, s*/
		pi, two_pi,	/* (2.0) (3.14...)			*/
		scale, 		/*  scale factor  */
		Tp = 1.0,	/* pulse period				*/
		Nc = 1.0,	/* number of pulse cycles		*/
		T = 20.0,	/* duration				*/
		sr = 100.0,	/* sample rate			 	*/
 		n  = 2.0,	/* rise-time value for pulse		*/
		phi = 0.0,	/* phase shift angle, radians		*/
		tau, 		/* decay time constant for the pulse	*/
		gain = 1.0,	/* update gain				*/
		*accel, *veloc, *displ, /* the pulse motion		*/
		*vc, *ac,	/* baseline correction records		*/
		GT, 		/* envelope normalization factor	*/
		velT, dspT, dspT_old; /* displacement and velocity at end */

	int	iter, maxIter = 50,	/* for baseline correction	*/
		p, p0 = 0,	/* number of initial points with no motion */
		points = 0,	/* number of points in the output file	*/
		sfrv = 0;	/*  scanf return value			*/

	void	exit();

	if ( argc == 1 ) {	// print usage
		fprintf(stderr," usage: \n");
		fprintf(stderr," pulse sr Tp Nc T pulse.dat \n");
        }

	if ( argc < 5 ) {	// interactive  input
		printf(" Sample rate:                 ");sfrv=scanf("%f", &sr );
		printf(" Pulse period (sec):          ");sfrv=scanf("%f", &Tp );
		printf(" Number of pulse cycles:      ");sfrv=scanf("%f", &Nc );
		printf(" Duration:                    ");sfrv=scanf("%f", &T );
		printf(" Output pulse record filename:");sfrv=scanf("%s", pulse_fn);
	} else {	// command-line input
		sr = atof(argv[1]);
		Tp = atof(argv[2]);
		Nc = atof(argv[3]);
		T  = atof(argv[4]);
		strcpy(pulse_fn , argv[5]);
		printf(" Tp=%f Nc=%f sr=%f T=%f fn=%s\n",Tp,Nc,sr,T,pulse_fn);
	}

	pi     = acos(-1.0);
	two_pi = 2.0*pi;

	points = floor(T*sr);
	p0     = floor(t0*sr);

 	tau = Nc*Tp/(2.0*n);	// decay time constant of pulse

	pulse_time =  2.0 * n * Nc * Tp;
	printf(" Pulse Time %.2f sec\n", pulse_time );

	if (1.0/(Tp*sr) > 0.051 )	// > 20 points / cycle 
	 printf("  warning: sample rate is too low for hydraulic actuators!\n");

	if (Tp < 15.0/sr) {
	 printf(" warning : sr=%f, Tp=%f, Tp must be greater than 15/sr\n",
				sr, Tp);
	 exit(1);
	}

	if ( 6.0*tau*n+t0 >= T ) {
	 printf(" warning: Tmax=%f, Nc*Tp/(2n)=%f, Tmax should be greater than 4 Np*Tp + t0\n", T, tau);
	 exit(1);
	}

	if (( pulse_fp = fopen(pulse_fn, "w" )) == NULL ) {
		printf ("  error: cannot open sine-sweep file '%s'\n",pulse_fn);
		exit(1);
	}
	
	// set phase of pulse such that terminal velocity is practically zero
	phi = -pi/2.0 + (n+1.0)*atan(tau*2*pi/Tp);

	// printf(" tau = %f phi = %f \n", tau, phi);	// debug

	accel  = vector(1,points);
	veloc  = vector(1,points);
	displ  = vector(1,points);
	vc     = vector(1,points);
	ac     = vector(1,points);

	for (p=1;p<=points;p++)	accel[p] = 0.0;
	for (p=1;p<=points;p++)	veloc[p] = 0.0;
	for (p=1;p<=points;p++)	displ[p] = 0.0;

	for (p=p0;p<=points;p++) {
		t = p/sr;
 		accel[p] = pow((t-t0)/tau,n) *
				exp(-(t-t0)/tau) * cos(two_pi*(t-t0)/Tp - phi);
	}

	GT = tau * 100;  /// fix this  /// 

// debug
//	printf("n=%f, T=%f tau=%f tg=%f gi=%e  GT=%e\n",
//			n, T, tau, exp(ln_gamma(n+1.0)), 1.0-gammainc(n+1.0,T), GT ); 
	// baseline correction for zero terminal velocity and displacement 
	for (p=1;p<=points;p++) {	// baseline correction functions
		t = p/sr;
		vc[p] = pow(((t-t0)/tau),n) * exp(-(t-t0)/tau) / GT;
		ac[p] = (exp(-(t-t0)/tau) * pow(((t-t0)/tau),(n-1.0)) * 
				(n*tau-t+t0)) / (tau*tau*GT); 
		if ( p < p0 ) vc[p] = ac[p] = 0.0;
	}

	velT = trapz(accel, 1.0/sr, points, 1 );
	dspT = trapz(accel, 1.0/sr, points, 2 );

	gain = 1.0;
	for (iter=1; iter<=maxIter; iter++) {	// baseline correction
		for (p=1; p<=points; p++) {
			accel[p] -= gain*( velT*vc[p] + dspT*ac[p] );
			if ( p < p0 ) accel[p] = 0.0;
		}
		dspT_old = dspT;
		velT = trapz(accel, 1.0/sr, points, 1 );
		dspT = trapz(accel, 1.0/sr, points, 2 );
        	if ( fabs(dspT) < 1e-5*Tp )	break;
		if ( fabs(dspT) > fabs(dspT_old) )	gain = gain/2.0;
	}

	veloc = cumTrapz(accel,1.0/sr,points);	// trapezoidal rule for veloc.

	displ = cumTrapz(veloc,1.0/sr,points);	// trapezoidal rule for displ.

	for (p=points-sr; p<=points; p++)	// end taper for last second
		displ[p] *= ((float)(points-p))/sr;

	scale = 2040.0 / maxAbsV(displ, points, &p);

	fprintf(pulse_fp, "%% pulse displacement data file \n");
	fprintf(pulse_fp, "%% pulse period     = %7.2f seconds \n", Tp );
	fprintf(pulse_fp, "%% number of cycles = %7.2f         \n", Nc );
	fprintf(pulse_fp, "%% duration         = %7.2f seconds \n", T );
	fprintf(pulse_fp, "%% sample rate      = %7.2f samples/second \n", sr );
	fprintf(pulse_fp, "%% pulse time       = %7.2f seconds \n",pulse_time);
	fprintf(pulse_fp, "%% iter = %d  dspT = %e  (baseline correction) \n", iter, dspT );

	for (p=1; p<=points; p++) 
		fprintf(pulse_fp, "%5.0f\n", displ[p]*scale );
//		fprintf(pulse_fp, "%15.5f %15.5f %15.5f\n",
//			accel[p]*scale, veloc[p]*scale, displ[p]*scale );
	
	printf(" %d points at %.2f samples per second in %.2f seconds\n",
					points, sr, points/sr );
	fclose (pulse_fp);

	free_vector(accel,1,points);
	free_vector(displ,1,points);
	free_vector(vc,1,points);
	free_vector(ac,1,points);
}
