/***************************************************************************
 Program scale.c - scales a time series recorded by an A/D
 converter and saved as integers between ADMIN and ADMAX. 
 Divides by the sensor sensitivity and deals with channel-to-channel skew.
 Optionaly corrects clipped data, smoothes data, and detrends data.  

 Clip Correction Types:
          0:  none
          3:  cubic polynomial
          5:  fifth order polynomial

 Detrending Types: 
    0:   none    
    1:   debias  
    2:   detrend 
    3:   baseline
    4:   first_pt
    5:   peak_peak

 To compile: make scale
 ... or ...
gcc -O -o scale scale.c NRutil.c  HPGmatrix.c HPGsignal.c HPGutil.c -lm

 To run: scale 'sensitivity file' 'unscaled file' 'scaled file' 'stats file'

 H.P. Gavin, Department of Civil Engineering, Duke University,   2024-02-11
****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "scale.h"
#include "../../HPGnumlib/HPGutil.h"
#include "../../HPGnumlib/NRutil.h"
#include "../../HPGnumlib/HPGmatrix.h"
#include "../../HPGnumlib/HPGsignal.h"

FILE    *fp_raw,    // pointer to the raw data file    
        *fp_scl,    // pointer to the scaled data file  
        *fp_stt;    // pointer to the statistics data file  


int main ( int argc, char *argv[] ) {

  char    ch,                   // dummy character      
          line[MAXL],           // a header line      
          title[MAXL],          // title of the sensitivity data file  
          xLabel[MAXL],         // label of X-axis data      
          yLabel[MAXL],         // label of Y-axis data      
          chnlLabel[16][MAXL],  // label for each channel  
          units[16][MAXL];      // sensitivity units of each channel  

  float   **r, **s,             // the raw and scaled data    
          sr=100.0, cr=1000.0,  // sample rate and conversion rate
	  t,                    // sample rate and current time    
          cf[24],               // correction factor for asynch samples  
          sensi[24],            // transducer sensitivity in volts/unit 
          min[24],max[24],      // minimum and maximum values    
          rms[24],              // root mean squared values    
          avg[24],              // average value of the time history  
          i_min=0,i_max=0,      // min and max of integrated signal  
          i_avg=0,i_rms=0,      // avg and rms of integrated signal  
          d_min=0,d_max=0,      // min and max of differentiated signal  
          d_avg=0,d_rms=0,      // avg and rms of differentiated signal  
          c[24][K+1],           // smoothing coefficients     
          s_f,                  // a smoothed and scaled data point  
          S[24],                // smoothing level between 0 and 1   
          wlo, whi,             // band-pass frequency limits 0..PI  
          range[16],            // voltage range of each channel        
          intg=0,               // an     integrated signal    
          diff=0,               // a         differentiated signal    
          UL, LL;               // upper and lower bounds of raw data  

  double  **N,                  // basis matrix for clip correction  
          Y[8];                 // data vector for clip correction  

  int     nScan=0,              // number of channel scans 
          nChnl,                // the number of channels    
          firstChnl,            // the first channel in the scan  
          lastChnl,             // the last channelin the scan    
          chn,                  // a channel number       
          len,                  // the length of a header line     
          C[24],                // clip correction poly order, 0, 3, or 5 
          D[24],                // the detrending type        
          integChnl,            // the channel to be integrated     
          diffrChnl,            // the channel to be differentiated   
          ok=1, pd=1, i, j, k, n, p;



  if (argc != 5) {
    color(1); color(31);
    printf("  usage:\n");
    color(41); color(33);
    printf("  scale <sensitivity file> <raw data file>");
    printf(" <scaled data file> <data stats file> ");
    printf("  argc = %d\n", argc);
    color(0); color(36);
    printf("\n");
     exit(1);
  }

  read_header ( argv, &nScan, &sr, &cr, &firstChnl, &lastChnl, &nChnl, range );

  read_sensitivity ( argv, title, sensi, C, D, S, 
                     firstChnl, lastChnl, xLabel, yLabel, chnlLabel, units,
                     &integChnl, &diffrChnl );

  if ( nChnl < 0 ) {
    printf ("%% file: %s ... \n", argv[2]);
    printf ("%% nScan: %8d  firstChnl: %2d  lastChnl: %2d   nChnl: %2d\n",
             nScan, firstChnl, lastChnl, nChnl);
    exit(1);
  }

  rewind (fp_raw);

// Allocate  Memory  For  The  Data  Matrices

  r = matrix( firstChnl, lastChnl, 1, nScan );
  s = matrix( firstChnl, lastChnl, 1, nScan );
  N = dmatrix(1,7,1,7);  // basis for clip correction  


  for ( chn=firstChnl; chn <= lastChnl; chn++ )
    cf[chn] = (float) (chn-firstChnl)*sr/BURST_RATE;

  for (i=1; i<=HEAD_LINES-1; i++) {  /* scan through the header */
    len = getLine ( fp_raw, MAXL, line );
    for (j=0; j<len-1; j++)
      putc ( line[j], fp_scl );
    putc('\n', fp_scl);
  }

// write sensitivity-related header metadata 
  fprintf(fp_scl, "%% sensor sensitivity file name: %s \n", argv[1] );
  fprintf(fp_scl, "%% %s \n", title );
  fprintf(fp_scl, "%% sensitivity:");
  for (chn=firstChnl; chn <= lastChnl; chn++) 
    fprintf(fp_scl,"%15.8e  ", sensi[chn]);
  fprintf(fp_scl,"\n");
  fprintf(fp_scl, "%% clip-corr. :");
  for (chn=firstChnl; chn <= lastChnl; chn++) 
    fprintf(fp_scl," %2d              ", C[chn]);
  fprintf(fp_scl,"\n");
  fprintf(fp_scl, "%% detrending :");
  for (chn=firstChnl; chn <= lastChnl; chn++) {
    switch ( D[chn] ) {
      case 0: fprintf(fp_scl," none            "); break;
      case 1: fprintf(fp_scl," debias          "); break;
      case 2: fprintf(fp_scl," detrend         "); break;
      case 3: fprintf(fp_scl," baseline        "); break;
      case 4: fprintf(fp_scl," firstPoint      "); break;
      case 5: fprintf(fp_scl," peak_peak       "); break;
    }
  }
  fprintf(fp_scl,"\n");
  fprintf(fp_scl, "%% smoothing  :");
  for (chn=firstChnl; chn <= lastChnl; chn++) 
    fprintf(fp_scl," %5.3f           ", S[chn]);
  fprintf(fp_scl,"\n");

  fprintf(fp_scl,"%%  time       ");
  for (chn=firstChnl; chn <= lastChnl; chn++)
    fprintf(fp_scl,"  chn %2d         ", chn );
  if ( integChnl >= firstChnl && integChnl <= lastChnl )
    fprintf(fp_scl," intg chn %2d ", integChnl );
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl )
    fprintf(fp_scl," diff chn %2d ", diffrChnl );
  fprintf(fp_scl,"\n");
  fprintf(fp_scl,"%%  seconds  ");
  for (chn=firstChnl; chn <= lastChnl; chn++)
    fprintf(fp_scl,"    %-11s  ",  units[chn] );
  fprintf(fp_scl,"\n");



// read in the raw integer data
  for ( p=1; p <= nScan; p++ ) {
    for ( chn=firstChnl; chn <= lastChnl; chn++ ) {
      ch=fscanf(fp_raw,"%f",&r[chn][p]);
      if ( (p>1) && (r[chn][p] == 0) ) r[chn][p] = r[chn][p-1]; // skipped conversion
    }
    while (( ch = getc(fp_raw)) != '\n') ;  // read to EOL 
  }

// correct clipped data using polynomial fits 
  for ( chn=firstChnl; chn <= lastChnl; chn++ ) { // loop over channels
    if ( C[chn] == 3 || C[chn] == 5 ) { // 3rd or 5th order polynomial 
      UL = LL = r[chn][1];              // upper and lower limits of data 
      for ( p=2; p <= nScan; p++ ) {
        if ( r[chn][p] > UL ) UL = r[chn][p];  // update UL 
        if ( r[chn][p] < LL ) LL = r[chn][p];  // update LL
      }
//    printf("| LL[%d] = %5.0f |  UL[%d] = %5.0f | \n",  chn,LL,chn,UL );

      p = 10; // no clipping correction for first 10 or last 10 nScan
      while ( p < nScan-10 ) {      // loop over all nScan channel scans
        p = p+1;

        if ( ( r[chn][p] == UL && r[chn][p+1] == UL ) || 
           ( r[chn][p] == LL && r[chn][p+1] == LL ) ) { // check for clip
           n = 0;
           while ( r[chn][p+n] == r[chn][p] )  ++n;
           if ( n > 2 ) {    // the data is clipped 
//           printf("| data clipped at point %6d | \n",  p );
       
             if (C[chn] == 3) {  // cubic polynomial clip correction 
               N[1][1] =   1.0;
               N[1][2] =  -2.0;
               N[1][3] =   4.0;
               N[1][4] =  -8.0;
               N[2][1] =   1.0;
               N[2][2] =  -1.0;
               N[2][3] =   1.0;
               N[2][4] =  -1.0;
               N[3][1] =   1.0;
               N[3][2] =  (n+1.0);
               N[3][3] =  (n+1.0)*(n+1.0);
               N[3][4] =  (n+1.0)*(n+1.0)*(n+1.0);
               N[4][1] =  1.0;
               N[4][2] =  (n+2.0);
               N[4][3] =  (n+2.0)*(n+2.0);
               N[4][4] =  (n+2.0)*(n+2.0)*(n+1.0);
         
               Y[1]    =  r[chn][p-2]; 
               Y[2]    =  r[chn][p-1]; 
               Y[3]    =  r[chn][p+n+1]; 
               Y[4]    =  r[chn][p+n+2]; 
         
               lu_dcmp ( N, 4, Y, 1, 1, &pd );  // solve for coef's 
       
               for (j=0; j<n; j++) 
                r[chn][p+j] = Y[1] + j*Y[2] + j*j*Y[3] + j*j*j*Y[4];
       
             } else {  // default fifth order polynomial clip correction 
       
               N[1][1] =    1.0;
               N[1][2] =   -3.0;
               N[1][3] =    9.0;
               N[1][4] =  -27.0;
               N[1][5] =   81.0;
               N[1][6] = -243.0;
               N[2][1] =    1.0;
               N[2][2] =   -2.0;
               N[2][3] =    4.0;
               N[2][4] =   -8.0;
               N[2][5] =   16.0;
               N[2][6] =  -32.0;
               N[3][1] =    1.0;
               N[3][2] =   -1.0;
               N[3][3] =    1.0;
               N[3][4] =   -1.0;
               N[3][5] =    1.0;
               N[3][6] =   -1.0;
               N[4][1] =    1.0;
               N[4][2] = (n+1.0);
               N[4][3] = (n+1.0)*(n+1.0);
               N[4][4] = (n+1.0)*(n+1.0)*(n+1.0);
               N[4][5] = (n+1.0)*(n+1.0)*(n+1.0)*(n+1.0);
               N[4][6] = (n+1.0)*(n+1.0)*(n+1.0)*(n+1.0)*(n+1.0);
               N[5][1] =    1.0;
               N[5][2] = (n+2.0);
               N[5][3] = (n+2.0)*(n+2.0);
               N[5][4] = (n+2.0)*(n+2.0)*(n+2.0);
               N[5][5] = (n+2.0)*(n+2.0)*(n+2.0)*(n+2.0);
               N[5][6] = (n+2.0)*(n+2.0)*(n+2.0)*(n+2.0)*(n+2.0);
               N[6][1] =    1.0;
               N[6][2] = (n+3.0);
               N[6][3] = (n+3.0)*(n+3.0);
               N[6][4] = (n+3.0)*(n+3.0)*(n+3.0);
               N[6][5] = (n+3.0)*(n+3.0)*(n+3.0)*(n+3.0);
               N[6][6] = (n+3.0)*(n+3.0)*(n+3.0)*(n+3.0)*(n+3.0);
               
               Y[1]    = r[chn][p-3]; 
               Y[2]    = r[chn][p-2]; 
               Y[3]    = r[chn][p-1]; 
               Y[4]    = r[chn][p+n+1]; 
               Y[5]    = r[chn][p+n+2]; 
               Y[6]    = r[chn][p+n+3]; 
       
               lu_dcmp ( N, 6, Y, 1, 1, &pd );  // solve for coef's 
    
               for (j=0; j<n; j++) 
                r[chn][p+j] = Y[1] + j*Y[2] + j*j*Y[3] +
                 j*j*j*Y[4] + j*j*j*j*Y[5] + j*j*j*j*j*Y[6];
       
             }      // end polynomial order switch
           }      // data is clipped
           p = p + n;
        }        // end clip check 
      }        // end data nScan loop
    }        // end C == 3 or C === 5 switch
  }        // end channel loop 

// scale data and correct for channel-to-channel skew 
  for ( p=1; p <= nScan; p++ ) {
    for ( chn=firstChnl; chn <= lastChnl; chn++ ) {

      r[chn][p] *= range[chn];    // voltage sensitivity 

      r[chn][p] /= DIGITAL_RANGE; // A/D sensitivity 
          
      s[chn][p] = r[chn][p];
      if ( p > 1 )        /* chnl-to-chnl skew    */
        s[chn][p] -= cf[chn]*(r[chn][p] - r[chn][p-1]);

      s[chn][p] /= sensi[chn];    /* sensor sensitivity   */
    }
  }

// detrend the data, as indicated 
  for (chn=firstChnl; chn <= lastChnl; chn++) {
    switch ( D[chn] ) {
      case 1: deBias     (s[chn], nScan); break;            // average value 
      case 2: deTrend    (s[chn], nScan); break;            // lst sqr line 
      case 3: baseLine   (s[chn], nScan,1,nScan,10); break; // base-line 
      case 4: firstPoint (s[chn], nScan, 10); break;        // first point 
      case 5: peakPeak   (s[chn], nScan); break;            // -min = +max 
    }
  }

  for ( chn=firstChnl; chn <= lastChnl; chn++ )
    for ( p=1; p <= nScan; p++ )  r[chn][p] = s[chn][p];

      /* set up symmetric low-pass filter coefficients */
  for (chn=firstChnl; chn <= lastChnl; chn++) {
    wlo = 0.0;
    whi = 0.15 + (PI-0.15) * ( 1.0 - S[chn] ); 
    c[chn][0] = ( whi - wlo ) / PI ;
    for (k=1; k <= K; k++)
      c[chn][k] = ( sin(whi*k) - sin(wlo*k) ) / (PI*k);
  }

                  /* filter the data */
  for ( p=1; p <= nScan; p++ ) {
    for ( chn=firstChnl; chn <= lastChnl; chn++ ) {
      if ( S[chn] > 0 && K < p && p < nScan-K ) {
        s_f = c[chn][0]*r[chn][p];  
        for (k=1; k < K; k++) s_f += c[chn][k]*( r[chn][p-k] + r[chn][p+k] );
        s[chn][p] = s_f;
      } else {
          s[chn][p] = r[chn][p];
      }
    }
  }

            /* find rms min & max values */
  for (chn=firstChnl; chn <= lastChnl; chn++) {
    min[chn] = max[chn] = rms[chn] = avg[chn] = 0.0;
    for ( p=1; p <= nScan; p++ ) {
      if ( p == 1 ) min[chn] = max[chn] = s[chn][p];
      if ( s[chn][p] < min[chn] ) min[chn] = s[chn][p];
      if ( s[chn][p] > max[chn] ) max[chn] = s[chn][p];
      rms[chn] += s[chn][p] * s[chn][p];
    }
    rms[chn] = sqrt ( rms[chn] / nScan );
  }

              /* save scaled data */
  for ( p=1; p <= nScan; p++ ) {
    fprintf(fp_scl,"%12.4e", (float) p/sr );
    for (chn=firstChnl; chn <= lastChnl; chn++)
      fprintf(fp_scl," %16.8e", s[chn][p] );

    /* integrate channel "integChnl" using trapezoidal rule */
    if ( integChnl >= firstChnl && integChnl <= lastChnl ) {
      chn = integChnl;
      if (p==1) {
        intg = 0.0;
        i_min = i_max = intg;
        i_rms = 0.0;
      } else {
        intg += ( s[chn][p-1] + s[chn][p] ) / (2.*sr);
        if ( intg > i_max )  i_max = intg;
        if ( intg < i_min )  i_min = intg;
        i_rms += intg*intg;
      }
      fprintf(fp_scl," %12.4e", intg );
    }

    /* differentiate channel "diffrChnl" using central diff. */
    if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) {
      chn = diffrChnl;
      if ( K+1 < p && p < nScan-K-1 )  /* smoothed */
        diff = sr*( s[chn][p+1] - s[chn][p-1] ) / 2.0;
      else         /* not smoothed  */
        diff = 0.0;
      if (p==1) {
        d_min = d_max = diff;
        d_rms = 0.0;
      } else {
        if ( diff > d_max )  d_max = diff;
        if ( diff < d_min )  d_min = diff;
        d_rms += diff*diff;
      }
        fprintf(fp_scl," %12.4e", diff );
    }

    fprintf(fp_scl,"\n");

  }

  if ( integChnl >= firstChnl && integChnl <= lastChnl ) 
    i_rms = sqrt ( i_rms / (float) nScan );
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) 
    d_rms = sqrt ( d_rms / (float) nScan );

  fprintf(fp_scl, "%% MINIMUM   ");
  for (chn=firstChnl; chn <= lastChnl; chn++)
    fprintf(fp_scl,"  %15.8e", min[chn] );
  if ( integChnl >= firstChnl && integChnl <= lastChnl ) 
    fprintf(fp_scl,"  %15.8e", i_min );
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) 
    fprintf(fp_scl,"  %15.8e", d_min );
  fprintf(fp_scl,"\n");
  fprintf(fp_scl, "%% MAXIMUM   ");
  for (chn=firstChnl; chn <= lastChnl; chn++)
    fprintf(fp_scl,"  %15.8e", max[chn] );
  if ( integChnl >= firstChnl && integChnl <= lastChnl ) 
    fprintf(fp_scl,"  %15.8e", i_max );
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) 
    fprintf(fp_scl,"  %15.8e", d_max );
  fprintf(fp_scl,"\n");
  fprintf(fp_scl, "%% R.M.S.    ");
  for (chn=firstChnl; chn <= lastChnl; chn++)
    fprintf(fp_scl,"  %15.8e", rms[chn] );
  if ( integChnl >= firstChnl && integChnl <= lastChnl ) 
    fprintf(fp_scl,"  %15.8e", i_rms );
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) 
    fprintf(fp_scl,"  %15.8e", d_rms );
  fprintf(fp_scl,"\n");

        /* print the data statistics to screen */
  color(1); color(33);
  printf("========================================");
  printf("======================================\n");
  color(1); color(32);
  printf(" channel             label ");
  printf("           min          max         rms \n");
  for (chn=firstChnl; chn <= lastChnl; chn++) {
    color(1); color(32);
    printf(" %2d ", chn );
    color(35);
    printf(" %26s ", chnlLabel[chn] );
    color(1); color(36);
    printf("    %10.3e   %10.3e  %10.3e \n", min[chn], max[chn], rms[chn] );
  }
  if ( integChnl >= firstChnl && integChnl <= lastChnl ) {
    color(1); color(32);
    printf(" %2d ", integChnl );
    color(35);
//  printf(" %12s ", chnlLabel[integChnl] );
    printf(" %-26s ", "integrated" );
    color(1); color(36);
    printf("    %10.3e   %10.3e  %10.3e \n", i_min, i_max, i_rms );
  }
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) {
    color(1); color(32);
    printf(" %2d ", diffrChnl );
    color(35);
//  printf(" %12s ", chnlLabel[diffrChnl] );
    printf(" %-26s ", "differentiated" );
    color(1); color(36);
    printf("    %10.3e   %10.3e  %10.3e \n", d_min, d_max, d_rms );
  }

        /* print the data statistics to file */
  fprintf(fp_stt,"%% data file: '%s'    %d  channel scans \n", argv[3], nScan );
  fprintf(fp_stt,"%% CHANNEL     MINIMUM      MAXIMUM     RMS \n");
  for (chn=firstChnl; chn <= lastChnl; chn++) {
    fprintf(fp_stt,"  %7d   %10.3e   %10.3e  %10.3e \n", 
          chn, min[chn], max[chn], rms[chn] );
  }
  if ( integChnl >= firstChnl && integChnl <= lastChnl ) 
    fprintf(fp_stt,"  intg %2d   %10.3e   %10.3e  %10.3e \n", 
          integChnl, i_min, i_max, i_rms );
  if ( diffrChnl >= firstChnl && diffrChnl <= lastChnl ) 
    fprintf(fp_stt,"  diff %2d   %10.3e   %10.3e  %10.3e \n", 
          diffrChnl, d_min, d_max, d_rms );
  fprintf(fp_stt,"%%_______________________________________________\n");

  color(0); color(36);

  fclose (fp_raw);
  fclose (fp_scl);
  fclose (fp_stt);

  free_matrix(r,firstChnl,lastChnl,1,nScan);
  free_matrix(s,firstChnl,lastChnl,1,nScan);
  free_dmatrix(N,1,7,1,7);

  return(1);
}


/*-----------------------------------------------------------------------------
READ_HEADER -  read data file header, open input/output files    7mar22
------------------------------------------------------------------------------*/
void read_header ( char *argv[], int *nScan, float *sr, float *cr, int *firstChnl, int *lastChnl, int *nChnl, float *range) 
{
  int  i, chn, ok;
  char  str[MAXL], ch;

  if (( fp_raw = fopen(argv[2], "r")) == NULL ) {
    color(1); color(31);
    printf("  error: cannot open raw data file '%s'\n", argv[2]);
    printf("  usage:\n");
    color(41); color(33);
    printf("  scale <sensitivity file> <raw data file>");
    printf(" <scaled data file> <data stats file> ");
    color(0); color(36);
    printf("\n");
    exit(1);
  }
  if (( fp_scl = fopen(argv[3], "w")) == NULL ) {
    color(1); color(31); 
    printf("  error: cannot open scaled data file '%s'\n", argv[3]);
    printf("  usage:\n");
    color(41); color(33);
    printf("  scale <sensitivity file> <raw data file>");
    printf(" <scaled data file> <data stats file> ");
    color(0); color(36);
    printf("\n");
    exit(1);
  }
  if (( fp_stt = fopen(argv[4], "a")) == NULL ) {
    color(1); color(31); 
    printf("  error: cannot open data stats file '%s'\n", argv[4]);
    printf("  usage:\n");
    color(41); color(33);
    printf("  scale <sensitivity file> <raw data file>");
    printf(" <scaled data file> <data stats file> ");
    color(0); color(36);
    printf("\n");
    exit(1);
  }


  // skip three lines
  for (i=1;i<=3;i++) while (( ch = getc(fp_raw)) != '\n') ;

// ComputerBoards DAC
// % 250 scans of channels 0 to 7 at 50 scans per second in 5.000 seconds
//ok=fscanf(fp_raw,"%s %d %s %s %s %d %s %d %s %f",
//                str, nScan, str,str,str, &firstChnl, str, &lastChnl, str, sr );
//  *nChnl = *firstChnl - *lastChnl + 1;

// HPADDA DAC
// % 1200 scans of 1 channels at 100 scans per second in 12.000 seconds
//ok=fscanf(fp_raw,"%s %d %s %s %d %s %s %f",
//                str, nScan, str, str, nChnl, str, str, sr );

// 2024-02-11 ...
// % 1200 scans of 2 channels at 100 sps and 2000 cps in 12.000 seconds
  ok=fscanf(fp_raw,"%s %d %s %s %d %s %s %f",
                  str, nScan, str, str, nChnl, str, str, sr, str, str, cr );
  
  *firstChnl = 0;  *lastChnl = *firstChnl + *nChnl - 1;

//fprintf(stdout,"\n\n nScan: %4d   firstnChnl:  %2d  lastChnl: %2d  nChnl:  %2d  sample rate: %.1f \n\n",
//            *nScan,  *firstChnl, *lastChnl, *nChnl, *sr); //*range ); // debug
  
  // skip three lines
  for (i=1;i<=3;i++) while (( ch = getc(fp_raw)) != '\n') ;
  ok=fscanf(fp_raw,"%s",str);   // comment character

  for (chn = 0; chn <= *nChnl; chn++) {
     ok=fscanf(fp_raw,"%f", &range[chn] );   // voltage measurement range
//   printf(" range[%d] = %f \n", chn, range[chn] ); 
  }

  rewind(fp_raw);
  return;
}


/*-----------------------------------------------------------------------------
READ_SENSITIVITY  - read sensitivity file 

    ---- SENSITIVITY  DATA  FILE  FORMAT ----  

Descriptive Title (one line only)
X label : "time, sec"
Y label : "volts"
integChnl :
diffrChnl :
Channel  Label             Sensitivity  V/Unit    DeClip  Detrend  Smooth  
===============================================================================
 0      "table accel."       0.004048   "cm/s/s"        0  3  0.1
 1      "table displ."       1.312      "cm"            0  3  0.1
 2      "platfm displ. 1"    1.000      "cm"            5  4  0.1
 3      "platfm displ. 2"    1.000      "cm"            5  4  0.1
 4      "platfm displ. 3"    1.000      "cm"            5  4  0.1
 5      "platfm accel. 1"    0.004048   "cm/s/s"        3  4  0.1
 6      "platfm accel. 2"    0.004048   "cm/s/s"        3  4  0.1
 7      "platfm accel. 3"    0.004048   "cm/s/s"        3  4  0.1

------------------------------------------------------------------------------*/
void read_sensitivity ( char *argv[], char *title,
      float *sensi, int *C, int *D, float *S,
      int firstChnl, int lastChnl,
      char *xLabel, char *yLabel,
      char (*chnlLabel)[MAXL], char (*units)[MAXL],
      int *integChnl, int *diffrChnl 
){
  FILE  *fp_sns;  /* pointer to the sensitivity data file  */
        int     i, chn, ok;
  char  str[MAXL];

  if ( ((fp_sns=fopen( argv[1],"r" )) == NULL )) {
   color(1); color(41); color(33);
   fprintf(stderr,"  error:  cannot open sensitivity data file '%s'   ",
            argv[1] );
   color(0); color(1); color(32);
   fprintf(stderr,"\n");
   fprintf(stderr,"\t\t---- SENSITIVITY  DATA  FILE  FORMAT ----\n");
   color(0); color(1); color(33);
   fprintf(stderr,"  Descriptive Title (one line only) \n");
   fprintf(stderr,"  X label : \"time, sec\"\n");
   fprintf(stderr,"  Y label : \"volts\"\n");
   fprintf(stderr," integChnl : -1 \n");
   fprintf(stderr," diffrChnl : -1 \n");
   fprintf(stderr,"  Channel       Label           Sensitivity  Units    DeClip Detrend Smooth\n");
   fprintf(stderr,"  =========================================================================\n");
   fprintf(stderr,"    :    \"descriptive label\"       :      \"units\" [0,3,5] [0-5]  [0-1] \n");
   color(0); color(36);
   exit(1);
  }

  for (i=0; i <= 15; i++)  sensi[i] = 1.0;  /* init. sensitivities */

  (void) getLine ( fp_sns, MAXL, title );
  scanLabel ( fp_sns, MAXL, xLabel, '"' );
  scanLabel ( fp_sns, MAXL, yLabel, '"' );
  scanLine ( fp_sns, MAXL, str, ':' );
  ok=fscanf ( fp_sns, "%d", integChnl );
  scanLine ( fp_sns, MAXL, str, ':' );
  ok=fscanf ( fp_sns, "%d", diffrChnl );
  (void) getLine ( fp_sns, MAXL, str );
  (void) getLine ( fp_sns, MAXL, str );
  (void) getLine ( fp_sns, MAXL, str );

  for (i=firstChnl; i <= lastChnl; i++) {
    ok=fscanf ( fp_sns,  "%d", &chn );
    scanLabel ( fp_sns, MAXL, chnlLabel[chn], '"' );
    ok=fscanf ( fp_sns, "%f", &sensi[chn] );
    scanLabel ( fp_sns, MAXL, units[chn], '"' );
    ok=fscanf ( fp_sns, "%d", &C[chn] );
    ok=fscanf ( fp_sns, "%d", &D[chn] );
    ok=fscanf ( fp_sns, "%f", &S[chn] );

    if ( chn < firstChnl || chn > lastChnl ) {
      color(1); color(41); color(33);
      fprintf(stderr," ... channel number %d ", chn );
      fprintf(stderr,"in sensitivity file '%s' is ", argv[1] );
      fprintf(stderr,"out of range ... ");
      color(0); color(1); color(33);
      fprintf(stderr,"\n\n");
      fprintf(stderr,"     start channel = %2d  \n", firstChnl );
      fprintf(stderr,"     stop channel  = %2d  \n", lastChnl );
      color(0); color(36);
      exit(0);
    }
    if ( S[chn] < 0 || S[chn] > 1 ) {
      color(1); color(41); color(33);
      fprintf(stderr," ... Smoothing level %4.2f ", S[chn] );
      fprintf(stderr," specified for channel %2d  \n", chn ); 
      fprintf(stderr,"     Smoothing level must be between 0 and 1 ");
      color(0); color(1); color(33);
      fprintf(stderr,"\n\n");
      fprintf(stderr,"     0 = no smoothing   \n");
      fprintf(stderr,"     1 = full smoothing \n");
      color(0); color(36);
      exit(0);
    }
    if ( D[chn] < 0 || D[chn] > 5 ) {
      color(1); color(41); color(33);
      fprintf(stderr," ... Detrending type %2d ", D[chn] );
      fprintf(stderr," specified for channel %2d  \n", chn ); 
      fprintf(stderr,"     Detrending type must be an integer 0 .. 5 ");
      color(0); color(1); color(33);
      fprintf(stderr,"\n\n");
      fprintf(stderr,"     0 = no detrending \n");
      fprintf(stderr,"     1 = de-bias     : subtract average value \n");
      fprintf(stderr,"     2 = de-trend    : subtract straight-line fit\n");
      fprintf(stderr,"     3 = baseline    : subract line through first and last points \n");
      fprintf(stderr,"     4 = first point : subract first point \n");
      fprintf(stderr,"     5 = peak-peak   : make max value = -min value \n");
      color(0); color(36);
      exit(0);
    }
  }

  color(1); color(37);
  printf("\n %s scaled via %s   \n", argv[2], argv[1] );
  color(1); color(33);
  printf("========================================");
  printf("======================================\n");
  color(1); color(32);
        printf(" %7s %17s  %16s  %-9s %6s %-8s %-6s\n",
    "channel", "label", "sensitivity", "units", "declip", "detrend", "smooth");
  for (chn= firstChnl; chn <= lastChnl; chn++) {
    color(32);
    printf(" %2d ", chn );
    color(35);
    printf(" %26s ", chnlLabel[chn] );
    color(36);
    printf("%12.3e  %-11s %1d  ", sensi[chn], units[chn], C[chn]);
    switch ( D[chn] ) {
      case 0: printf("  none     "); break;
      case 1: printf("  debias   "); break;
      case 2: printf("  detrend  "); break;
      case 3: printf("  baseline "); break;
      case 4: printf("  first_pt "); break;
      case 5: printf("  peak_peak"); break;
            }
    printf(" %5.2f\n", S[chn]);
  }

  fclose(fp_sns);

  return;
}
