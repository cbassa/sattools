#./bin/bash

# Loop into every directory
for obsdir in *; do 
	if [ -d $obsdir ]; then
		echo $obsdir
		cd $obsdir
		# check if observation is calibrated
		if [ ! -e "out.dat" ]; then
			CALFRAME=1
			while true; do
				# calibrate Nth image captured
				st_calibrate.sh $CALFRAME
				# check if no calibration data was output
				if [ ! -e "out.dat" ]; then
					 read -p "Delete observation?" yn
					 case $yn in
						  [Yy]* ) cd ..; rm -rf $obsdir; break;;
						  [Nn]* ) CALFRAME=$(($CALFRAME+1)); continue;;
						  * ) echo "Please answer y or n.";;
					 esac
				else
					cd ..
					break
				fi
			done
		else
			cd ..
		fi
	fi
done
