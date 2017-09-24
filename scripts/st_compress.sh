#!/bin/bash

source ~/.bashrc

# Settings
PGMDIR=/dev/shm
N=250
# If no images are found during this period (seconds), relaunch capture process
WAIT=20
COUNT=0
# Start automatically capture process or not
AUTOSTART=0

# Default camera to start imaging if no schedule is available
CAMERADEV="/dev/video0"
CAPTUREPID=0

# Check obsdir exists
if [ ! -d $ST_OBSDIR ]; then
	mkdir $ST_OBSDIR
fi
if [ ! -d $ST_OBSDIR/control ]; then
	mkdir $ST_OBSDIR/control
fi


# if autostart force a restart
if [ $AUTOSTART == 1 ]; then
	echo "restart" >$ST_OBSDIR/control/state.txt
	STATE="restart"
else
	echo "stop" >$ST_OBSDIR/control/state.txt
	STATE="stop"
fi

export CAMERADEV=`cat $ST_OBSDIR/control/camera.txt | awk '{print $((7))}'`

echo "Status: "$STATE

# For ever loop
while true; do
	
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
		#			 cp $ST_OBSDIR/control/scale.txt .
		cp $ST_OBSDIR/control/camera.txt .
		export CAMERADEV=`cat $ST_OBSDIR/control/camera.txt | awk '{print $((7))}'`
		# Remove old captured frames
		echo "Removing all captured frames"
		#		ls -1 $PGMDIR/img*.pgm | awk '{printf("sudo rm -rf %s\n",$1)}' | sh
		#ls -1 $PGMDIR/img*.pgm | awk '{printf("rm -rf %s\n",$1)}' | sh
		find $PGMDIR -type f -name './img*.pgm' | awk '{printf("rm -rf %s\n",$1)}' | sh
	fi

	# kill capture process just when scheduler sends stop signal
	if [ $STATE == "stop" ]; then
		if [ $(($CAPTUREPID)) != 0 ]; then
			echo ""
			echo "Stopping capture process"
			kill -9 $CAPTUREPID
			CAPTUREPID=0
		fi
	fi

	# Get number of captured frames
	NFILES=`ls -1 $PGMDIR/img*.pgm 2>/dev/null | wc -l`
	#NFILES=`find $PGMDIR -type f -name 'img*.pgm' |wc -l`
	# If enough, process
	if [ $NFILES -ge $N ]; then
		COUNT=0
		echo ""
		echo "Compressing $N captured frames"

		# Start point
		M=`ls -1 $PGMDIR/img*.pgm | head -n1 | sed -e "s/[^0-9]*//g"` 

		# Run pgm2fits
#			pgm2fits -p $PGMDIR/img -w 720 -h 576 -s $M -n $N 
		pgm2fits -p $PGMDIR/img -w 720 -h 576 -s $M -n $N >/dev/null

		# Run viewer
		viewer `ls -1 2*.fits | tail -n1`
		cp avg.pgm $ST_OBSDIR

		# Remove files
		echo ""
		echo "Removing $N captured frames"
#		ls -1 $PGMDIR/img*.pgm | head -n$N | awk '{printf("sudo rm -rf %s\n",$1)}' | sh
		ls -1 $PGMDIR/img*.pgm | head -n$N | awk '{printf("rm -rf %s\n",$1)}' | sh
		#find $PGMDIR -type f -name 'img*.pgm' | head -n$N | awk '{printf("rm -rf %s\n",$1)}' | sh

#		echo "Finished"
	else
	  # There are not enough captured frames
		# Launch capture process if state is not stop
		# if time passes with still no images re-launch capture process
    if [ $STATE != "stop" ]; then
			echo ""
    	echo "Waiting for frames. Status: "$STATE
			COUNT=$(($COUNT+1))
			if [ $COUNT -ge $WAIT ];	then
				COUNT=0
				echo ""
				echo "No frames found, restarting capture script"
#				sh $ST_DATADIR/scripts/st_capture.sh /dev/video-$CAMERA &
				sh $ST_DATADIR/scripts/st_capture.sh $CAMERADEV &
				sleep 1
				CAPTUREPID=`pgrep -o -x ffmpeg`
			fi
			# if restarting then relaunch capture script now
			if [ $STATE == "restart" ]; then
				if [ $(($CAPTUREPID)) == 0 ]; then
					echo "Restarting capture script"
#					sh $ST_DATADIR/scripts/st_capture.sh /dev/video-$CAMERA &
					sh $ST_DATADIR/scripts/st_capture.sh $CAMERADEV &
					sleep 1
					CAPTUREPID=`pgrep -o -x ffmpeg`
					COUNT=0
				fi
			fi
    else
		# we are stopped, check for bogus capture process
#    	echo "Status: "$STATE
			if [ $(($CAPTUREPID)) != 0 ]; then
				echo "Stopping capture process"
				kill -9 $CAPTUREPID
				CAPTUREPID=0
			fi
    fi
	fi
  # Sleep
  sleep 1
done
