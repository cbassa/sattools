#!/bin/bash


# Check obsdir exists
if [ ! -d $ST_OBSDIR ]; then
	mkdir $ST_OBSDIR
fi
if [ ! -d $ST_OBSDIR/control ]; then
	mkdir $ST_OBSDIR/control
fi

# look for observations in every subdir
# catenate each one to lastobs.txt in obsdir
find $ST_OBSDIR -iname 'observations.txt' | awk '{printf ("cat %s\n",$1)}' | sh > $ST_OBSDIR/lastobs.txt
#  and rename each original obs files to *.txt.used
find $ST_OBSDIR -iname 'observations.txt' | awk '{printf ("mv %s %s.used\n",$1,$1)}' | sh


# catenate lastobs.txt to allobs.txt and count the new observations
cat $ST_OBSDIR/lastobs.txt >> $ST_OBSDIR/allobs.txt
# count total observations
NOBS=`cat $ST_OBSDIR/lastobs.txt | wc -l`
NALLOBS=`cat $ST_OBSDIR/allobs.txt | wc -l`

echo "$NOBS new observations in lastobs.txt, total observations in allobs.txt: $NALLOBS"

