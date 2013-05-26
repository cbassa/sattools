# Makefile: http://www.eng.hawaii.edu/Tutor/Make/

# Compiling flags
CFLAGS = -O3 -Wno-unused-result

# Linking flags
LFLAGS = -lm -lcpgplot -lpgplot -lX11 -fno-backslash -lpng -lqfits -lwcs_c -lgsl -lgslcblas -ljpeg -lexif

# Compilers
CC = gcc
F77 = gfortran

all:
	make satfit uk2iod rde2iod viewer residuals tleinfo satmap satorbit runsched fitskey fitsheader satid skymap addwcs reduce wcsfit plotfits pgm2fits

devel:	
	make faketle imgstat

faketle: faketle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o faketle faketle.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

imgstat: imgstat.o
	$(CC) -o imgstat imgstat.o $(LFLAGS)


satfit: satfit.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o
	$(F77) -o satfit satfit.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o $(LFLAGS)

uk2iod: uk2iod.o
	$(CC) -o uk2iod uk2iod.o $(LFLAGS)

rde2iod: rde2iod.o
	$(CC) -o rde2iod rde2iod.o $(LFLAGS)

viewer: viewer.o
	$(CC) -o viewer viewer.o $(LFLAGS)

residuals: residuals.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o residuals residuals.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

tleinfo: tleinfo.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o tleinfo tleinfo.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

satmap: satmap.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o satmap satmap.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

satorbit: satorbit.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o satorbit satorbit.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

runsched: runsched.o
	$(CC) -o runsched runsched.o $(LFLAGS)

fitskey: fitskey.o
	$(CC) -o fitskey fitskey.o -lqfits

fitsheader: fitsheader.o
	$(CC) -o fitsheader fitsheader.o -lqfits

satid: satid.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o satid satid.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

skymap: skymap.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o skymap skymap.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

reduce: reduce.o
	$(F77) -o reduce reduce.o $(LFLAGS)

addwcs: addwcs.o
	$(F77) -o addwcs addwcs.o $(LFLAGS)

wcsfit: wcsfit.o
	$(F77) -o wcsfit wcsfit.o $(LFLAGS)

plotfits: plotfits.o
	$(F77) -o plotfits plotfits.o $(LFLAGS)

pgm2fits: pgm2fits.o
	$(F77) -o pgm2fits pgm2fits.o $(LFLAGS) 

clean:
	rm -f *.o
	rm -f *~
