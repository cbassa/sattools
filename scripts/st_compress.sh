#!/bin/bash

# Settings
PGMDIR=/dev/shm
N=250
# If no images are found during thisperiod, relaunch capture process
WAIT=20
COUNT=0
# Start automatically capture process or not
AUTOSTART=0

# Default camera to start imaging if no schedule is available
CAMERA="S25H"

# if autostart force a restart
if [ $AUTOSTART == 1 ]; then
	echo "restart" >$ST_OBSDIR/control/state.txt
	STATE="restart"
else
	echo "stop" >$ST_OBSDIR/control/state.txt
	STATE="stop"
fi

export CAMERA=`cat $ST_OBSDIR/control/camera.txt | awk '{print $((1))}'`

# For ever loop
while true; do
	# Get number of files
	NFILES=`ls -1 $PGMDIR/img*.pgm 2>/dev/null | wc -l`
	
	# Get state
	export STATE=`cat $ST_OBSDIR/control/state.txt`
	# If enough, process
	if [ $NFILES -ge $N ]; then
		 
		# Create new directory
		if [ $STATE == "restart" ]; then
			 export DIR=`date -u +%FT%T | sed -e "s/-//g" -e "s/\://g" -e "s|T|/|g"`
			 mkdir -p $ST_OBSDIR/$DIR
			 cd $ST_OBSDIR/$DIR
			 echo "Moving to $ST_OBSDIR/$DIR"
			 echo "observing" >$ST_OBSDIR/control/state.txt
			 cp $ST_OBSDIR/control/position.txt .
#			 cp $ST_OBSDIR/control/scale.txt .
			 cp $ST_OBSDIR/control/camera.txt .
       export CAMERA=`cat $ST_OBSDIR/control/camera.txt | awk '{print $((1))}'`
		fi

		# Process only if not stopped
		if [ $STATE != "stop" ]; then
			echo "Compressing captured frames"

			# Start point
			M=`ls -1 $PGMDIR/img*.pgm | head -n1 | sed -e "s/[^0-9]*//g"` 

			# Run pgm2fits
#			pgm2fits -p $PGMDIR/img -w 720 -h 576 -s $M -n $N 
			pgm2fits -p $PGMDIR/img -w 720 -h 576 -s $M -n $N >/dev/null

			# Run viewer
			viewer `ls -1 2*.fits | tail -n1`
			cp avg.pgm $ST_OBSDIR
		else
			kill -9 $CAPTUREPID
			CAPTUREPID=0
		fi

		# Remove files
		echo "Removing captured frames"
#		ls -1 $PGMDIR/img*.pgm | head -n$N | awk '{printf("sudo rm -rf %s\n",$1)}' | sh
		ls -1 $PGMDIR/img*.pgm | head -n$N | awk '{printf("rm -rf %s\n",$1)}' | sh

		echo "Finished"
	else
	  # There are not enough images available
		# Launch capture process
		# if time passes with still no images re-launch capture process
    if [ $STATE != "stop" ]; then
    	echo "Waiting for images. Status: "$STATE
			COUNT=$(($COUNT+1))
			if [ $COUNT -ge $WAIT ];	then
				COUNT=0
				echo "No images found, restarting capture script"
				sh $ST_DATADIR/scripts/st_capture.sh /dev/video-$CAMERA &
				sleep 1
				CAPTUREPID=`pgrep -o -x ffmpeg`
			fi
			if [ $STATE == "restart" ]; then
				if [ $(($CAPTUREPID)) == 0 ]; then
					sh $ST_DATADIR/scripts/st_capture.sh /dev/video-$CAMERA &
					sleep 1
					CAPTUREPID=`pgrep -o -x ffmpeg`
				fi
			fi
    else
    	echo "Status: "$STATE
			if [ $(($CAPTUREPID)) != 0 ]; then
				kill -9 $CAPTUREPID
				CAPTUREPID=0
			fi
    fi
	fi
  # Sleep
  sleep 1
done
