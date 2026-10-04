# SGP4/SDP4 library (git submodule) and sattools' additions to it
SGP4DIR = ../extern/dundee_sgp4
SGP4LIB = $(SGP4DIR)/libsgdp4.a
SGP4 = sattools_util.o $(SGP4LIB)

# Compiling flags
CFLAGS = -O3 -Wno-unused-result -I$(SGP4DIR)

# Linking flags
LFLAGS = -lm -lcpgplot -lpgplot -lX11 -lwcs -lgsl -lgslcblas -lpng 

# Compilers
CC = gcc
F77 = gfortran

all: check-sgp4
	make addwcs angular calibrate dec2sex faketle fitsheader fitskey imgstat jpg2fits jpgstack measure pgm2fits plotfits pstrack rde2iod reduce residuals runsched satfit satid satmap satorbit sex2dec skymap tle2ole tleinfo uk2iod stviewer wcsfit deproject slewto waitfor pass detect launchtle propagate fakeiod csv2tle normal posmatch posvel xyz2tle mvtle confirm allnight

lite: check-sgp4
	make faketle residuals satfit satmap satorbit tleinfo uk2iod pass launchtle propagate skymap allnight

geolon: geolon.o $(SGP4)
	$(CC) -o ../bin/geolon geolon.o $(SGP4) -lm

allnight: allnight.o
	$(CC) -o ../bin/allnight allnight.o -lm

selectiod: selectiod.o
	$(CC) -o ../bin/selectiod selectiod.o -lm

planscan: planscan.o $(SGP4)
	$(CC) -o ../bin/planscan planscan.o $(SGP4) -lm

tle2rv: tle2rv.o $(SGP4)
	$(CC) -o ../bin/tle2rv tle2rv.o $(SGP4) -lm

mvtle: mvtle.o $(SGP4)
	$(CC) -o ../bin/mvtle mvtle.o $(SGP4) -lm

rv2tle: rv2tle.o $(SGP4)
	$(CC) -o ../bin/rv2tle rv2tle.o $(SGP4) -lm

xyz2tle: xyz2tle.o versafit.o dsmin.o simplex.o $(SGP4)
	$(CC) -o ../bin/xyz2tle xyz2tle.o versafit.o dsmin.o simplex.o $(SGP4) -lm

csv2tle: csv2tle.o $(SGP4)
	$(CC) -o ../bin/csv2tle csv2tle.o $(SGP4) -lm

fakeiod: fakeiod.o $(SGP4)
	$(CC) -o ../bin/fakeiod fakeiod.o $(SGP4) -lm

posvel: posvel.o $(SGP4)
	$(CC) -o ../bin/posvel posvel.o $(SGP4) -lm

normal: normal.o $(SGP4)
	$(CC) -o ../bin/normal normal.o $(SGP4) -lm

vadd: vadd.o $(SGP4)
	$(CC) -o ../bin/vadd vadd.o $(SGP4) -lm

posmatch: posmatch.o $(SGP4)
	$(CC) -o ../bin/posmatch posmatch.o $(SGP4) -lm

propagate: propagate.o $(SGP4)
	$(CC) -o ../bin/propagate propagate.o $(SGP4) -lm

detect: detect.o forward.o reverse.o
	$(F77) -o ../bin/detect detect.o forward.o reverse.o -lm $(LFLAGS) -lqfits

confirm: confirm.o forward.o reverse.o
	$(F77) -o ../bin/confirm confirm.o forward.o reverse.o -lm $(LFLAGS) -lqfits

autodetect: autodetect.o
	$(F77) -o ../bin/autodetect autodetect.o -lm $(LFLAGS)

launchtle: launchtle.o $(SGP4)
	$(CC) -o ../bin/launchtle launchtle.o $(SGP4) -lm

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

pstrack: pstrack.o forward.o reverse.o $(SGP4)
	$(F77) -o ../bin/pstrack pstrack.o forward.o reverse.o $(SGP4) $(LFLAGS) -lqfits

faketle: faketle.o $(SGP4)
	$(CC) -o ../bin/faketle faketle.o $(SGP4) -lm

imgstat: imgstat.o
	$(CC) -o ../bin/imgstat imgstat.o -lm -lqfits

satfit: satfit.o versafit.o dsmin.o simplex.o $(SGP4)
	$(F77) -o ../bin/satfit satfit.o versafit.o dsmin.o simplex.o $(SGP4) $(LFLAGS)

uk2iod: uk2iod.o
	$(CC) -o ../bin/uk2iod uk2iod.o -lm

rde2iod: rde2iod.o
	$(CC) -o ../bin/rde2iod rde2iod.o -lm

stviewer: stviewer.o
	$(CC) -o ../bin/stviewer stviewer.o -lm -lqfits

residuals: residuals.o $(SGP4)
	$(CC) -o ../bin/residuals residuals.o $(SGP4) -lm -lwcs

tleinfo: tleinfo.o $(SGP4)
	$(CC) -o ../bin/tleinfo tleinfo.o $(SGP4) -lm

satmap: satmap.o $(SGP4)
	$(F77) -o ../bin/satmap satmap.o $(SGP4) $(LFLAGS)

satorbit: satorbit.o $(SGP4)
	$(F77) -o ../bin/satorbit satorbit.o $(SGP4) $(LFLAGS)

runsched: runsched.o
	$(CC) -o ../bin/runsched runsched.o -lm

tle2ole: tle2ole.o
	$(CC) -o ../bin/tle2ole tle2ole.o -lm

fitskey: fitskey.o
	$(CC) -o ../bin/fitskey fitskey.o -lqfits

fitsheader: fitsheader.o
	$(CC) -o ../bin/fitsheader fitsheader.o -lqfits

satid: satid.o forward.o reverse.o $(SGP4)
	$(F77) -o ../bin/satid satid.o forward.o reverse.o $(SGP4) $(LFLAGS) -lqfits

skymap: skymap.o $(SGP4)
	$(F77) -o ../bin/skymap skymap.o $(SGP4) $(LFLAGS)

pass: pass.o $(SGP4)
	$(CC) -o ../bin/pass pass.o $(SGP4) -lm

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
	-$(MAKE) -C $(SGP4DIR) clean

# SGP4/SDP4 library and element set readers, from the dundee_sgp4 submodule
$(SGP4LIB): FORCE
	$(MAKE) -C $(SGP4DIR) lib CC="$(CC)"

check-sgp4:
	@test -f $(SGP4DIR)/sgdp4h.h || { echo "extern/dundee_sgp4 is empty: run 'git submodule update --init'"; exit 1; }

FORCE:

.PHONY: all lite clean check-sgp4 FORCE
