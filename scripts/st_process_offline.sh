#!/bin/bash

# Loop over files
for file in 2*.fits; do
    # Run sextractor
    sextractor $file -c $ST_DATADIR/sextractor/default.sex 2>/dev/null
    mv test.cat $file.cat

    # Run addwcs
    addwcs -f $file -r test.fits

    # Run satid
    satid $file $file.png/png 2>/dev/null 

    # Move calibrated file
    mv $file.cat $file.cal $file.id $file $file.png png/
done
