#!/bin/bash

if [ -z $1 ]; then
    N=10
else
		echo "Using image $1 to calibrate"
    N=$1
		rm test.fits
fi

if [ ! -e "test.fits" ]; then
    ls -1 2*.fits | head -n$N | tail -n1 | awk '{printf("cp %s test.fits\n",$1)}' | sh
fi
sextractor test.fits -c $ST_DATADIR/sextractor/default.sex
mv test.cat test.fits.cat

mkdir cal png

plotfits test.fits
wcsfit
