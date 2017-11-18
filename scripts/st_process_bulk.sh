#./bin/bash

# Loop into every directory
for obsdir in *; do 
	if [ -d $obsdir ]; then
		echo $obsdir
		cd $obsdir
		# check if observation is not already processed
		N=`ls -1 ./2*.fits 2>/dev/null | wc -l`
		if [ $N -ge 1 ]; then
			# process observation frames in directory
			st_process_offline.sh $1
		fi
		cd ..
	fi
done


