CC= gcc-mp-4.5
FC= gfortran
CLIBS= -lm
CFLAGS= -O -g -Wall
FFLAGS= -O -fno-automatic
DEFINES=-DNO_JPEG -DNO_LCMS

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) -c $<

tphot: tphot.o nitfread.o rwfits.o psf2d.o fitlm.o linsolve.o sortc.o
	$(CC) $(CFLAGS) tphot.o nitfread.o rwfits.o \
	 psf2d.o fitlm.o linsolve.o sortc.o -lm -o tphot
