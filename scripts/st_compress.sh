#!/bin/bash

# Settings
PGMDIR=/dev/shm
N=250

# Force a restart
echo "restart" >$ST_OBSDIR/control/state.txt

# For ever loop
while true; do
    # Get number of files
    NFILES=`ls -1 $PGMDIR/img*.pgm 2>/dev/null | wc -l`

    # If enough, process
    if [ $NFILES -ge $N ]; then
	# Get state
	export STATE=`cat $ST_OBSDIR/control/state.txt`
    
	# Create new directory
	if [ $STATE == "restart" ]; then
	    export DIR=`date -u +%FT%T | sed -e "s/-//g" -e "s/\://g" -e "s|T|/|g"`
	    mkdir -p $ST_OBSDIR/$DIR
	    cd $ST_OBSDIR/$DIR
	    echo "Moving to $ST_OBSDIR/$DIR"
	    echo "observing" >$ST_OBSDIR/control/state.txt
	    cp $ST_OBSDIR/control/position.txt .
	    cp $ST_OBSDIR/control/camera.txt .
	fi

	# Start point
	M=`ls -1 $PGMDIR/img*.pgm | head -n1 | sed -e "s/[^0-9]*//g"` 

	# Run pgm2fits
	pgm2fits -p $PGMDIR/img -w 720 -h 576 -s $M -n $N 

	# Remove files
	ls -1 $PGMDIR/img*.pgm | head -n$N | awk '{printf("sudo rm -rf %s\n",$1)}' | sh

	# Run viewer
	viewer `ls -1 2*.fits | tail -n1`
	cp avg.pgm $ST_OBSDIR
	echo "Finished"
    fi

    # Sleep
    sleep 1
done
