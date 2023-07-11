# HPGdaac-xtra

Programs to support digital data acquisition using **HPGdaac**
High Performance Graphical data acquisition and control

---------------------------------

## Installation 


1. clone software from github e.g., to your ~/Code/ directory

```
mkdir ~/Code
cd ~/Code
git clone https://github.com/hpgavin/HPGnumlib 
git clone https://github.com/hpgavin/HPGdaac-xtra
```

2. make and make install

```
cd ~/Code/HPGdaac-xtra
make clean
make
sudo make install
```

---------------------------------

## Usage

**HPGdaac-xtra** contains a set of sixteen command line programs written in C to support digitial data acquisition from analog sensors using **HPGdaac**.
There are two kinds of programs:  

1. Programs to generate waveform data files to be converted to analog signals by **HPGdaac**
2. Programs to process  waveform data files measured and recorded by **HPGdaac**


### Programs to generate waveform data files to be converted to analog signals by **HPGdaac**

* **bwrand.c** - Generate a vector of gaussean band-limited random noise 
```
  input:  sr      : sample rate (samples per second) 
          offset  : offset from zero
          ampl    : ampliutde (percent full scale)
          f_lo    : low  frequency limit of band-limited noise
          f_hi    : high frequency limit of band-limited noise
          T       : time duration of the data
          Seed    : seed of the random number generator
          outFile : output data file

  output: u       : random output signal
```
usage
```
bwrand  sr offset  f_lo  f_hi  T   seed   outfile
```

* **chirpq.c** - Generate a file of square waveform data with changing frequency and amplitude
```
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

  output: u       : square wave data
```
usage
```
chirps sr offset ao af fo ff T p q outFile 
```

* **chirps.c** - Generate a file of sinusoidal waveform data with changing frequency and amplitude
```
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

  output: u       : sine wave data
```
usage
```
chirps sr offset ao af fo ff T p q outFile 
```

* **chirpt.c** - Generate a file of triangular waveform data with changing frequency and amplitude
```
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

  output: u       : triangle wave data
```
usage
```
chirpt sr offset ao af fo ff T p q outFile
```

* **ktrand.c** - synthetic earthquake ground motion with the Kannai Tajimi spectrum
```
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
```
usage
``` 
ktrand  sr RMS fg zg aa tt T seed datafile
```

* **pulse.c** - Generate a file with pulse-like waveform data.  

```
  input:  sr      : sample rate, samples per second 
          Tp      : period of the pulse
          Nc      : number of cycles in the pulse
          T       : duration of the record
          Vp      : pulse amplitude ( > 0ber generator
          file    : name of output file
```
usage

```
pulse sr Tp Nc T pulse.dat 

```

* **ramp.c** - Generate a file of ramp data with constant velocity 

usage
```
ramp
```

------------------------------


### Programs to process waveform data files measured and recorded by **HPGdaac**

* **baseline.c** - Perform baseline correction on waveform data
```
  input:  infile   : input data file
  output: output   : output data file
```
usage
```
baseline <infile> <outfile>
```

* **glue.c** -  glues data file A to data file B, column-wise.

Data files A and B are multi-column data files, which have a set of 'header' 
lines starting with a # or % character.  Data file A can have more rows than
data file B, but can not have less rows than data file B.   

The 'glued' file will use the header of data file A and will place the columns
of data file B after the last column of data file A.  If data file A is longer
than data file B, the additional rows of data will be zero's.  

... currently, file B must be a one-column data file ... 

usage
```
glue fileA fileB fileC
```

* **limits.c**  -  scans a data file and keeps track of max and min 
usage
```
limits <filename>
```

* **resample.c** -  Re-samples a multi-column time history data file.

usage
```
resample in_file  out_file
```

* **scale.c** - scales a time series recorded by **HPGdaac** and saved as integers between ADMIN and ADMAX. 
 Divides by the sensor sensitivity and deals with channel-to-channel skew.
 Optionaly corrects clipped data, smoothes data, and detrends data.  

```
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
```
usage
```
scale 'sensor configuration file' 'raw data file' 'scaled data file' 'stats file'
```
Sensor configuration file format
```      
Descriptive Title (one line only)
X label : "time, sec" 
Y label : "volts"
integChnl :
diffrChnl :
Channel  Label             Sensitivity  V/Unit    DeClip  Detrend  Smooth  
===============================================================================
 0  "table accel."       0.004048   "cm/s/s"  0  3  0.1
 1      "table displ."       1.312      "cm"            0  3  0.1
 2      "platfm displ. 1"    1.000      "cm"            5  4  0.1
 3      "platfm displ. 2"    1.000      "cm"            5  4  0.1
 4      "platfm displ. 3"    1.000      "cm"            5  4  0.1
 5      "platfm accel. 1"    0.004048   "cm/s/s"        3  4  0.1
 6      "platfm accel. 2"    0.004048   "cm/s/s"        3  4  0.1
 7      "platfm accel. 3"    0.004048   "cm/s/s"        3  4  0.1
```

* **skew.c**  -  calculate the channel to channel skew of an AD device

* **xfer.c** - calculates frequency response functions H1, H2, and Hv, and the
associated power spectra Gu, Gy, and Guy using data sets x & y. 
Averaging and Windowing are Incorporated. 

usage
```
xfer <data file>
```
------------------------------------
