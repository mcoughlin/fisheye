CC= cc
FC= gfortran
CLIBS= -lm
# CFLAGS= -O -g -Wall
CFLAGS= -O -g
FFLAGS= -O -fno-automatic
DCRAWDEF=-DNO_JPEG -DNO_LCMS

mdcraw.o: mdcraw.c
	$(CC) $(CFLAGS) -c mdcraw.c $(DCRAWDEF) $(LIBS) -o mdcraw.o

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) -c $<

linear: linear.o invert.o
	$(FC) $(FFLAGS) linear.o invert.o $(FLIBS) -o linear

median: median.o
	$(CC) $(CFLAGS) median.o -lm -o median

colmerge: colmerge.o
	$(CC) $(CFLAGS) colmerge.o -o colmerge

fitrot: fitrot.o minimrq.o invert_c.o
	$(CC) $(CFLAGS) fitrot.o minimrq.o invert_c.o -lm -o fitrot 

tphot: tphot.o nitfread.o rwfits.o psf2d.o fitlm.o linsolve.o sortc.o
	$(CC) $(CFLAGS) tphot.o nitfread.o rwfits.o \
	 psf2d.o fitlm.o linsolve.o sortc.o -lm -o tphot

pgm2fits: pgm2fits.o
	$(CC) $(CFLAGS) pgm2fits.o -lm -o pgm2fits

mdcraw: mdcraw.o
	$(CC) $(CFLAGS) mdcraw.o $(DCRAWDEF) -lm -o mdcraw

catstars: catstars.o
	$(CC) $(CFLAGS) catstars.o -o catstars

catsect: catsect.o
	$(CC) $(CFLAGS) catsect.o -o catsect

cam2sky: cam2sky.o rotten.o
	$(CC) cam2sky.o rotten.o -o cam2sky -lm

clean:
	rm *.o tphot mdcraw median linear fitshdr fitrot colmerge pgm2fits catstars catsect

all:
	make linear
	make median
	make colmerge
	make fitrot
	make tphot
	make pgm2fits
	make mdcraw
	make fitshdr
	make catstars
	make catsect
	make cam2sky
