#!/bin/bash

# Take into account that for this script to work as expected it should be run as ". st_storeobs.sh" with
# a dot at the beginning of the command line

OBSTIME=`pwd | awk -F"/" '{print $(NF)}'`
DATE=`pwd | awk -F"/" '{print $(NF-1)}'`

echo $OBSTIME
echo $DATE

cd .. | sh

if [ ! -d $ST_OBSDIR/stored ]; then
	mkdir $ST_OBSDIR/stored
fi

if [ ! -d $ST_OBSDIR/stored/$DATE ]; then
	mkdir $ST_OBSDIR/stored/$DATE
fi

cd ..
mv $OBSTIME $ST_OBSDIR/stored/$DATE

