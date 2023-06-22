SHELL = /bin/bash

DIR_C := ./src
DIR_N := ../HPGnumlib
DIR_O := ./obj

CC := gcc
CFLAGS := -O
LFLAGS := -lm 

all : baseline  bwrand  chirpq  chirps  chirpt   deskip  glue  ktrand  limits   pulse  ramp  resample  scale  skew  xfer 

$(DIR_O)/%.o : $(DIR_N)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(DIR_O)/%.o : $(DIR_C)/%.c
	$(CC) $(CFLAGS) -c $< -o $@


baseline : $(DIR_O)/baseline.o  $(DIR_O)/HPGutil.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)
	
bwrand : $(DIR_O)/bwrand.o  $(DIR_O)/HPGutil.o  $(DIR_O)/NRutil.o   $(DIR_O)/HPGsignal.o   $(DIR_O)/HPGmatrix.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)
	
chirpq : $(DIR_O)/chirpq.o  $(DIR_O)/HPGutil.o  $(DIR_O)/NRutil.o  
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

chirps : $(DIR_O)/chirps.o  $(DIR_O)/HPGutil.o  $(DIR_O)/NRutil.o  
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

chirpt : $(DIR_O)/chirpt.o  $(DIR_O)/HPGutil.o  $(DIR_O)/NRutil.o  
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

deskip : $(DIR_O)/deskip.o  $(DIR_O)/HPGutil.o  $(DIR_O)/NRutil.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

glue : $(DIR_O)/glue.o  $(DIR_O)/HPGutil.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

ktrand : $(DIR_O)/ktrand.o  $(DIR_O)/HPGutil.o  $(DIR_O)/NRutil.o  $(DIR_O)/HPGmatrix.o $(DIR_O)/HPGsignal.o $(DIR_O)/HPGrandom.o  $(DIR_O)/HPGcholesky.o  
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

limits : $(DIR_O)/limits.o  
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

pulse : $(DIR_O)/pulse.o $(DIR_O)/NRutil.o  $(DIR_O)/HPGsignal.o  $(DIR_O)/HPGmatrix.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

ramp : $(DIR_O)/ramp.o  
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

resample : $(DIR_O)/resample.o   $(DIR_O)/HPGutil.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

scale : $(DIR_O)/scale.o   $(DIR_O)/NRutil.o   $(DIR_O)/HPGutil.o  $(DIR_O)/HPGmatrix.o  $(DIR_O)/HPGsignal.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

skew : $(DIR_O)/skew.o   
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

xfer : $(DIR_O)/xfer.o  $(DIR_O)/NRutil.o   $(DIR_O)/HPGutil.o  $(DIR_O)/HPGmatrix.o  $(DIR_O)/NRcomplex.o  $(DIR_O)/HPGsignal.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)


install :
	mv baseline bwrand chirpq chirps chirpt deskip glue ktrand limits pulse ramp resample scale skew xfer  /usr/local/bin

clean :
	rm $(DIR_O)/*.o 

