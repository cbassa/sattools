#!/bin/bash

# For ever loop
while true; do
    # Get number of files
    NFILES=`ls -1 2*.fits 2>/dev/null | wc -l`

    # If enough, process
    if [ $NFILES -gt 0 ]; then
	
	# Get file name
	FILE=`ls -1 2*.fits | head -n1`
	sleep 1

	# Run sextractor
	sextractor $FILE -c $ST_DATADIR/sextractor/default.sex 2>/dev/null
	mv test.cat $FILE.cat

	# Run addwcs
	addwcs -f $FILE -r test.fits 

	# Run satid
	satid $FILE $FILE.png/png 2>/dev/null

	
	# Move calibrated file
	mv $FILE $FILE.cat $FILE.cal $FILE.png $FILE.id $FILE.det png/ 2>/dev/null
    fi
    sleep 1
done
