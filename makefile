# Makefile: http://www.eng.hawaii.edu/Tutor/Make/

# Compiling flags
CFLAGS = #-O3 -Wno-unused-result

# Linking flags
LFLAGS = -lm -lcpgplot -lpgplot -lX11 -lwcs -lgsl -lgslcblas -lpng

# Compilers
CC = gcc
F77 = gfortran

all:
	make addwcs angular calibrate dec2sex faketle fitsheader fitskey imgstat jpg2fits jpgstack measure pgm2fits plotfits pstrack rde2iod reduce residuals runsched satfit satid satmap satorbit sex2dec skymap tle2ole tleinfo uk2iod stviewer wcsfit deproject slewto waitfor pass detect launchtle propagate fakeiod csv2tle normal posmatch posvel xyz2tle mvtle confirm allnight

lite:
	make faketle residuals satfit satmap satorbit tleinfo uk2iod pass launchtle propagate skymap allnight

geolon: geolon.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o geolon geolon.o sgdp4.o satutl.o deep.o ferror.o -lm

allnight: allnight.o
	$(CC) -o allnight allnight.o -lm

selectiod: selectiod.o
	$(CC) -o selectiod selectiod.o -lm

planscan: planscan.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o planscan planscan.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)	

tle2rv: tle2rv.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o tle2rv tle2rv.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)	

mvtle: mvtle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o mvtle mvtle.o sgdp4.o satutl.o deep.o ferror.o -lm

rv2tle: rv2tle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o rv2tle rv2tle.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

xyz2tle: xyz2tle.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o
	$(CC) -o xyz2tle xyz2tle.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o -lm

csv2tle: csv2tle.o satutl.o ferror.o
	$(CC) -o csv2tle csv2tle.o satutl.o ferror.o -lm

fakeiod: fakeiod.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o fakeiod fakeiod.o sgdp4.o satutl.o deep.o ferror.o -lm

posvel: posvel.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o posvel posvel.o sgdp4.o satutl.o deep.o ferror.o -lm

normal: normal.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o normal normal.o sgdp4.o satutl.o deep.o ferror.o -lm

vadd: vadd.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o vadd vadd.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

posmatch: posmatch.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o posmatch posmatch.o sgdp4.o satutl.o deep.o ferror.o -lm

propagate: propagate.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o propagate propagate.o sgdp4.o satutl.o deep.o ferror.o -lm

detect: detect.o forward.o reverse.o
	$(F77) -o detect detect.o forward.o reverse.o -lm $(LFLAGS) -lqfits

confirm: confirm.o forward.o reverse.o
	$(F77) -o confirm confirm.o forward.o reverse.o -lm $(LFLAGS) -lqfits

autodetect: autodetect.o
	$(F77) -o autodetect autodetect.o -lm $(LFLAGS)

launchtle: launchtle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o launchtle launchtle.o sgdp4.o satutl.o deep.o ferror.o -lm

slewto: slewto.o
	$(CC) -o slewto slewto.o -lm

waitfor: waitfor.o
	$(CC) -o waitfor waitfor.o -lm

deproject: deproject.o forward.o reverse.o
	$(F77) -o deproject deproject.o forward.o reverse.o $(LFLAGS) -ljpeg -lqfits

jpgstack: jpgstack.o
	$(CC) -o jpgstack jpgstack.o -ljpeg

angular: angular.o forward.o reverse.o
	$(CC) -o angular angular.c forward.o reverse.o -lm -lwcs

dec2sex: dec2sex.o 
	$(CC) -o dec2sex dec2sex.c -lm

sex2dec: sex2dec.o 
	$(CC) -o sex2dec sex2dec.c -lm

calibrate: calibrate.o forward.o
	$(F77) -o calibrate calibrate.o forward.o $(LFLAGS) -lqfits

measure: measure.o reverse.o
	$(F77) -o measure measure.o reverse.o $(LFLAGS) -lqfits

jpg2fits: jpg2fits.o
	$(CC) -o jpg2fits jpg2fits.o -lm -lqfits -ljpeg

pstrack: pstrack.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o
	$(F77) -o pstrack pstrack.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o $(LFLAGS) -lqfits

faketle: faketle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o faketle faketle.o sgdp4.o satutl.o deep.o ferror.o -lm

imgstat: imgstat.o
	$(CC) -o imgstat imgstat.o -lm -lqfits

satfit: satfit.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o forward.o
	$(F77) -o satfit satfit.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o forward.o $(LFLAGS)

uk2iod: uk2iod.o
	$(CC) -o uk2iod uk2iod.o -lm

rde2iod: rde2iod.o
	$(CC) -o rde2iod rde2iod.o -lm

stviewer: stviewer.o
	$(CC) -o stviewer stviewer.o -lm -lqfits

residuals: residuals.o sgdp4.o satutl.o deep.o ferror.o forward.o
	$(CC) -o residuals residuals.o sgdp4.o satutl.o deep.o ferror.o forward.o -lm -lwcs

tleinfo: tleinfo.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o tleinfo tleinfo.o sgdp4.o satutl.o deep.o ferror.o -lm

satmap: satmap.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o satmap satmap.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

satorbit: satorbit.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o satorbit satorbit.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

runsched: runsched.o
	$(CC) -o runsched runsched.o -lm

fitskey: fitskey.o
	$(CC) -o fitskey fitskey.o -lqfits

fitsheader: fitsheader.o
	$(CC) -o fitsheader fitsheader.o -lqfits

satid: satid.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o
	$(F77) -o satid satid.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o $(LFLAGS) -lqfits

skymap: skymap.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o skymap skymap.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

pass: pass.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o pass pass.o sgdp4.o satutl.o deep.o ferror.o -lm

reduce: reduce.o forward.o reverse.o
	$(F77) -o reduce reduce.o forward.o reverse.o $(LFLAGS) -lqfits

addwcs: addwcs.o forward.o reverse.o
	$(F77) -o addwcs addwcs.o forward.o reverse.o $(LFLAGS) -lqfits

wcsfit: wcsfit.o forward.o reverse.o
	$(F77) -o wcsfit wcsfit.o forward.o reverse.o $(LFLAGS) -lqfits

plotfits: plotfits.o forward.o reverse.o
	$(F77) -o plotfits plotfits.o forward.o reverse.o $(LFLAGS) -lqfits

pgm2fits: pgm2fits.o
	$(F77) -o pgm2fits pgm2fits.o $(LFLAGS) -lqfits

clean:
	rm -f *.o
	rm -f *~
