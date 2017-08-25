#!/bin/bash

# Settings
PGMDIR=/dev/shm
N=250
WAIT=10
COUNT=0

# Force a restart
echo "restart" >$ST_OBSDIR/control/state.txt
STATE="restart"

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
			 cp $ST_OBSDIR/control/scale.txt .
			 cp $ST_OBSDIR/control/camera.txt .
       export CAMERA=`cat $ST_OBSDIR/control/camera.txt`
		fi

		# Process only if not stopped
		if [ $STATE != "stop" ]; then
			# Start point
			M=`ls -1 $PGMDIR/img*.pgm | head -n1 | sed -e "s/[^0-9]*//g"` 

			# Run pgm2fits
			pgm2fits -p $PGMDIR/img -w 720 -h 576 -s $M -n $N 

			# Run viewer
			viewer `ls -1 2*.fits | tail -n1`
			cp avg.pgm $ST_OBSDIR
		fi

		# Remove files
		ls -1 $PGMDIR/img*.pgm | head -n$N | awk '{printf("sudo rm -rf %s\n",$1)}' | sh

		echo "Finished"
  # if not enough images are still available wait to re-launch capture process
	else
	echo "Waiting for images"
    if [ $STATE == "observing" ]; then
			COUNT=$(($COUNT+1))
			if [ $COUNT -ge $WAIT ];	then
				COUNT=0
				echo "No images found, restarting capture script"
				sh $ST_DATADIR/scripts/st_capture.sh /dev/video-$CAMERA&
	#		else
	#			echo $COUNT
			fi
    fi
	fi
  # Sleep
  sleep 1
done
