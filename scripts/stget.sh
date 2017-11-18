#!/bin/bash

# Get date
DATE=`date +%Y%m%d_%H%M%S`

# Get cookie
wget --post-data='identity=yourlogin&password=yourpassword' --cookies=on --keep-session-cookies --save-cookies=/tmp/cookies.txt 'https://www.space-track.org/ajaxauth/login' -o /tmp/stget.log

# Get data
wget  --keep-session-cookies --load-cookies=/tmp/cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/EPOCH/%3Enow-30/format/3le' -O $ST_TLEDIR/catalog.tle
dos2unix $ST_TLEDIR/catalog.tle
sed -i -e "s/^1     /1 0000/g" -e "s/^2     /2 0000/g" -e "s/^1    /1 000/g" -e "s/^2    /2 000/g" -e "s/^1   /1 00/g" -e "s/^2   /2 00/g" -e "s/^1  /1 0/g" -e "s/^2  /2 0/g" $ST_TLEDIR/catalog.tle
cp $ST_TLEDIR/catalog.tle $ST_TLEDIR/${DATE}_catalog.txt

# Get classfd
wget http://www.prismnet.com/~mmccants/tles/classfd.zip --no-check-certificate -O $ST_TLEDIR/classfd.zip
unzip -o $ST_TLEDIR/classfd.zip
dos2unix $ST_TLEDIR/classfd.tle
cp $ST_TLEDIR/classfd.tle $ST_TLEDIR/${DATE}_classfd.txt
#mv $HOME/classfd.tle $ST_TLEDIR/classfd.tle
rm $ST_TLEDIR/classfd.zip

# Get inttles
wget http://www.prismnet.com/~mmccants/tles/inttles.zip --no-check-certificate -O $ST_TLEDIR/inttles.zip
unzip -o $ST_TLEDIR/inttles.zip
dos2unix $ST_TLEDIR/inttles.tle
cp $ST_TLEDIR/inttles.tle $ST_TLEDIR/${DATE}_inttles.txt
#mv $HOME/inttles.tle $ST_TLEDIR/inttles.tle
rm $ST_TLEDIR/inttles.zip

#rm $HOME/login

# Create bulk file
cat $ST_TLEDIR/classfd.tle $ST_TLEDIR/catalog.tle >$ST_TLEDIR/bulk.tle
