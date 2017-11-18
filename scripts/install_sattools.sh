#!/bin/bash

echo "Step 1.1: make directory"
sleep 1
mkdir -p $HOME/code/c/satellite/
cd $HOME/code/c/satellite

echo "Step 1.2: clone sattools repository (may take a while)"
sleep 1
#git clone https://github.com/cbassa/sattools.git
git clone https://github.com/fmederos/sattools.git

echo "Step 1.3: make"
sleep 1
cd $HOME/code/c/satellite/sattools
make

echo "Step 1.4: clone strf repository"
sleep 1
cd $HOME/code/c/satellite
git clone https://github.com/cbassa/strf.git

echo "Step 1.5: make"
sleep 1
cd $HOME/code/c/satellite/strf
make

echo "Step 1.6: download classfd.tle"
sleep 1
mkdir -p $HOME/code/c/satellite/sattools/tle
cd $HOME/code/c/satellite/sattools/tle
wget -c https://www.prismnet.com/~mmccants/tles/classfd.zip
unzip classfd.zip
dos2unix classfd.tle
rm classfd.zip

echo "Step 2.1: set environment variables"
sleep 1
cd
echo "export PGPLOT_DIR=/usr/local/src/pgplot-5.2.2" >>$HOME/.xsessionrc
echo "export ST_COSPAR=4171" >>$HOME/.xsessionrc
echo "export ST_DATADIR=$HOME/code/c/satellite/sattools" >>$HOME/.xsessionrc
echo "export ST_TLEDIR=$HOME/code/c/satellite/sattools/tle" >>$HOME/.xsessionrc
echo "export ST_OBSDIR=$HOME/satobs" >>$HOME/.xsessionrc
mkdir $ST_OBSDIR
mkdir $ST_OBSDIR/control

echo "Step 2.2: set path"
sleep 1
echo "PATH=$HOME/code/c/satellite/sattools:$HOME/code/c/satellite/sattools/scripts:$HOME/software/strf:\$PATH" >>$HOME/.profile

echo "Final step: run"
echo "source $HOME/.profile"
