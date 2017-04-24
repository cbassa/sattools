#!/bin/bash

if [ ! -e "test.fits" ]; then
    ls -1 2*.fits | head -n10 | tail -n1 | awk '{printf("cp %s test.fits\n",$1)}' | sh
fi
sextractor test.fits -c $ST_DATADIR/sextractor/default.sex
mv test.cat test.fits.cat

mkdir cal png

plotfits test.fits
wcsfit
