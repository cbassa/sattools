# Compiling flags
CFLAGS = -O3 -Wno-unused-result

# Linking flags
LFLAGS = -lm -lcpgplot -lpgplot -lX11 -lwcs -lgsl -lgslcblas -lpng -lgiza -L../qfits/lib

# Compilers
CC = gcc
F77 = gfortran

all:
	make addwcs angular calibrate dec2sex faketle fitsheader fitskey imgstat jpg2fits jpgstack measure pgm2fits plotfits pstrack rde2iod reduce residuals runsched satfit satid satmap satorbit sex2dec skymap tle2ole tleinfo uk2iod stviewer wcsfit deproject slewto waitfor pass detect launchtle propagate fakeiod csv2tle normal posmatch posvel xyz2tle mvtle confirm allnight

lite:
	make faketle residuals satfit satmap satorbit tleinfo uk2iod pass launchtle propagate skymap allnight

geolon: geolon.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/geolon geolon.o sgdp4.o satutl.o deep.o ferror.o -lm

allnight: allnight.o
	$(CC) -o ../bin/allnight allnight.o -lm

selectiod: selectiod.o
	$(CC) -o ../bin/selectiod selectiod.o -lm

planscan: planscan.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/planscan planscan.o sgdp4.o satutl.o deep.o ferror.o -lm

tle2rv: tle2rv.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/tle2rv tle2rv.o sgdp4.o satutl.o deep.o ferror.o -lm

mvtle: mvtle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/mvtle mvtle.o sgdp4.o satutl.o deep.o ferror.o -lm

rv2tle: rv2tle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/rv2tle rv2tle.o sgdp4.o satutl.o deep.o ferror.o -lm

xyz2tle: xyz2tle.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o
	$(CC) -o ../bin/xyz2tle xyz2tle.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o -lm

csv2tle: csv2tle.o satutl.o ferror.o
	$(CC) -o ../bin/csv2tle csv2tle.o satutl.o ferror.o -lm

fakeiod: fakeiod.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/fakeiod fakeiod.o sgdp4.o satutl.o deep.o ferror.o -lm

posvel: posvel.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/posvel posvel.o sgdp4.o satutl.o deep.o ferror.o -lm

normal: normal.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/normal normal.o sgdp4.o satutl.o deep.o ferror.o -lm

vadd: vadd.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/vadd vadd.o sgdp4.o satutl.o deep.o ferror.o -lm

posmatch: posmatch.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/posmatch posmatch.o sgdp4.o satutl.o deep.o ferror.o -lm

propagate: propagate.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/propagate propagate.o sgdp4.o satutl.o deep.o ferror.o -lm

detect: detect.o forward.o reverse.o
	$(F77) -o ../bin/detect detect.o forward.o reverse.o -lm $(LFLAGS) -lqfits

confirm: confirm.o forward.o reverse.o
	$(F77) -o ../bin/confirm confirm.o forward.o reverse.o -lm $(LFLAGS) -lqfits

autodetect: autodetect.o
	$(F77) -o ../bin/autodetect autodetect.o -lm $(LFLAGS)

launchtle: launchtle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/launchtle launchtle.o sgdp4.o satutl.o deep.o ferror.o -lm

slewto: slewto.o
	$(CC) -o ../bin/slewto slewto.o -lm

waitfor: waitfor.o
	$(CC) -o ../bin/waitfor waitfor.o -lm

deproject: deproject.o forward.o reverse.o
	$(F77) -o ../bin/deproject deproject.o forward.o reverse.o $(LFLAGS) -ljpeg -lqfits

jpgstack: jpgstack.o
	$(CC) -o ../bin/jpgstack jpgstack.o -ljpeg

angular: angular.o forward.o reverse.o
	$(CC) -o ../bin/angular angular.c forward.o reverse.o -lm -lwcs

dec2sex: dec2sex.o 
	$(CC) -o ../bin/dec2sex dec2sex.c -lm

sex2dec: sex2dec.o 
	$(CC) -o ../bin/sex2dec sex2dec.c -lm

calibrate: calibrate.o forward.o
	$(F77) -o ../bin/calibrate calibrate.o forward.o $(LFLAGS) -lqfits

measure: measure.o reverse.o
	$(F77) -o ../bin/measure measure.o reverse.o $(LFLAGS) -lqfits

jpg2fits: jpg2fits.o
	$(CC) -o ../bin/jpg2fits jpg2fits.o -lm -lqfits -ljpeg

pstrack: pstrack.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o
	$(F77) -o ../bin/pstrack pstrack.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o $(LFLAGS) -lqfits

faketle: faketle.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/faketle faketle.o sgdp4.o satutl.o deep.o ferror.o -lm

imgstat: imgstat.o
	$(CC) -o ../bin/imgstat imgstat.o -lm -lqfits

satfit: satfit.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o
	$(F77) -o ../bin/satfit satfit.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o $(LFLAGS)

uk2iod: uk2iod.o
	$(CC) -o ../bin/uk2iod uk2iod.o -lm

rde2iod: rde2iod.o
	$(CC) -o ../bin/rde2iod rde2iod.o -lm

stviewer: stviewer.o
	$(CC) -o ../bin/stviewer stviewer.o -lm -lqfits

residuals: residuals.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/residuals residuals.o sgdp4.o satutl.o deep.o ferror.o -lm -lwcs

tleinfo: tleinfo.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/tleinfo tleinfo.o sgdp4.o satutl.o deep.o ferror.o -lm

satmap: satmap.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o ../bin/satmap satmap.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

satorbit: satorbit.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o ../bin/satorbit satorbit.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

runsched: runsched.o
	$(CC) -o ../bin/runsched runsched.o -lm

tle2ole: tle2ole.o
	$(CC) -o ../bin/tle2ole tle2ole.o -lm

fitskey: fitskey.o
	$(CC) -o ../bin/fitskey fitskey.o -lqfits

fitsheader: fitsheader.o
	$(CC) -o ../bin/fitsheader fitsheader.o -lqfits

satid: satid.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o
	$(F77) -o ../bin/satid satid.o sgdp4.o satutl.o deep.o ferror.o forward.o reverse.o $(LFLAGS) -lqfits

skymap: skymap.o sgdp4.o satutl.o deep.o ferror.o
	$(F77) -o ../bin/skymap skymap.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

pass: pass.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o ../bin/pass pass.o sgdp4.o satutl.o deep.o ferror.o -lm

reduce: reduce.o forward.o reverse.o
	$(F77) -o ../bin/reduce reduce.o forward.o reverse.o $(LFLAGS) -lqfits

addwcs: addwcs.o forward.o reverse.o
	$(F77) -o ../bin/addwcs addwcs.o forward.o reverse.o $(LFLAGS) -lqfits

wcsfit: wcsfit.o forward.o reverse.o
	$(F77) -o ../bin/wcsfit wcsfit.o forward.o reverse.o $(LFLAGS) -lqfits

plotfits: plotfits.o forward.o reverse.o
	$(F77) -o ../bin/plotfits plotfits.o forward.o reverse.o $(LFLAGS) -lqfits

pgm2fits: pgm2fits.o
	$(F77) -o ../bin/pgm2fits pgm2fits.o $(LFLAGS) -lqfits

clean:
	rm -f *.o
	rm -f *~
	rm -f ../bin/*
