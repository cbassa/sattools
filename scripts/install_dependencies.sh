#!/bin/bash

echo "Step 1.1: install dependencies"
sleep 1
apg-get update
apt-get install  emacs gfortran libpng-dev libx11-dev libjpeg-dev libexif-dev git dos2unix sextractor

echo "Step 1.2: goto /usr/local/src"
sleep 1
cd /usr/local/src

echo "Step 2.1: download pgplot"
sleep 1
wget -c ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz

echo "Step 2.2: unpack pgplot"
sleep 1
gunzip -c pgplot5.2.tar.gz | tar xvf -

echo "Step 2.3: create pgplot directory"
sleep 1
mkdir -p /usr/local/src/pgplot-5.2.2

echo "Step 2.4: select drivers"
sleep 1
# Selecting PNDRIV, PSDRIV and XWDRIV
sed -e "s/! PNDRIV/  PNDRIV/g" -e "s/! PSDRIV/  PSDRIV/g" -e "s/! XWDRIV/  XWDRIV/g" pgplot/drivers.list >pgplot-5.2.2/drivers.list

echo "Step 2.5: create makefile"
sleep 1
cd /usr/local/src/pgplot-5.2.2
../pgplot/makemake ../pgplot linux g77_gcc

echo "Step 2.6: adjusting makefile"
sleep 1
sed -i -e "s/FCOMPL=g77/FCOMPL=gfortran/g" makefile
sed -i -e "s/FFLAGC=-u -Wall -fPIC -O/FFLAGC=-ffixed-form -ffixed-line-length-none -u -Wall -fPIC -O/g" makefile
sed -i -e "s|pndriv.o : ./png.h ./pngconf.h ./zlib.h ./zconf.h|pndriv.o : |g" makefile

echo "Step 2.7: run make"
sleep 1
make
make cpg

echo "Step 2.8: place libraries and header files"
sleep 1
rm -rf /usr/local/lib/libpgplot.a /usr/local/lib/libcpgplot.a /usr/local/lib/libpgplot.so /usr/local/include/cpgplot.h
ln -s /usr/local/src/pgplot-5.2.2/libpgplot.a /usr/local/lib/
ln -s /usr/local/src/pgplot-5.2.2/libpgplot.so /usr/local/lib/
ln -s /usr/local/src/pgplot-5.2.2/libcpgplot.a /usr/local/lib/
ln -s /usr/local/src/pgplot-5.2.2/cpgplot.h /usr/local/include/

echo "Step 2.9: clean up"
sleep 1
rm -rf /usr/local/src/pgplot5.2.tar.gz /usr/local/src/pgplot

echo "Step 3.1: download qfits"
sleep 1
cd /usr/local/src
wget -c ftp://ftp.eso.org/pub/qfits/qfits-5.2.0.tar.gz

echo "Step 3.2: unpack qfits"
sleep 1
gunzip -c qfits-5.2.0.tar.gz | tar xvf -

echo "Step 3.3: fix xmemory.c"
sleep 1
cd /usr/local/src/qfits-5.2.0
chmod +w src/xmemory.c
sed -i -e "s/swapfd = open(fname, O_RDWR | O_CREAT);/swapfd = open(fname, O_RDWR | O_CREAT, 0644);/g" src/xmemory.c

echo "Step 3.4: configure and make"
sleep 1
./configure
make
make install

echo "Step 3.5: clean up"
sleep 1
rm /usr/local/src/qfits-5.2.0.tar.gz

echo "Step 4.1: download wcslib-2.9"
sleep 1
cd /usr/local/src
wget -c http://www.epta.eu.org/~bassa/wcslib-2.9.tar

echo "Step 4.2: unpack wcslib"
sleep 1
tar -xvf wcslib-2.9.tar

echo "Step 4.3: compile wcslib"
sleep 1
cd /usr/local/src/wcslib-2.9/C/
make clean
rm libwcs_c.a
make

echo "Step 4.4: place libraries and header files"
sleep 1
rm -rf /usr/local/lib/libwcs_c.a /usr/local/include/proj.h /usr/local/include/cel.h
ln -s /usr/local/src/wcslib-2.9/C/libwcs_c.a /usr/local/lib/
ln -s /usr/local/src/wcslib-2.9/C/proj.h /usr/local/include/
ln -s /usr/local/src/wcslib-2.9/C/cel.h /usr/local/include/

echo "Step 4.5: clean up"
sleep 1
rm -rf /usr/local/src/wcslib-2.9.tar

echo "Step 5.1: download gsl"
sleep 1
cd /usr/local/src
wget -c ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz

echo "Step 5.2: unpack gsl"
sleep 1
gunzip -c gsl-1.16.tar.gz | tar xvf -

echo "Step 5.3: configure, make, make install"
sleep 1
cd /usr/local/src/gsl-1.16/
./configure
make 
make install

echo "Step 5.4: clean up"
sleep 1
rm -rf /usr/local/src/gsl-1.16.tar.gz

echo "Step 6.1: set ld.so.conf"
sleep 1
echo "include /etc/ld.so.conf.d/*.conf" >/etc/ld.so.conf
echo "/usr/local/lib" >>/etc/ld.so.conf
ldconfig

echo "Step 6.1: download fftw"
sleep 1
cd /usr/local/src
wget http://www.fftw.org/fftw-3.3.4.tar.gz

echo "Step 6.2: unpack fftw"
sleep 1
gunzip -c fftw-3.3.4.tar.gz | tar xvf -

echo "Step 6.3: configure, make,make install"
sleep 1
cd /usr/local/src/fftw-3.3.4
./configure --enable-float
make 
make install

echo "Step 6.4: clean up"
sleep 1
rm -rf /usr/local/src/fftw-3.3.4.tar.gz

echo "Step 7.1: download ffmpeg"
sleep 1
cd /usr/local/src
wget http://ffmpeg.org/releases/ffmpeg-snapshot.tar.bz2

echo "Step 7.2: unpack ffmpeg"
sleep 1
bzip2 -cd ffmpeg-snapshot.tar.bz2 | tar xvf -

echo "Step 7.3: patch pgmenc.c"
sleep 1
cd /usr/local/src/ffmpeg/libavcodec
wget https://dl.dropboxusercontent.com/u/52579487/pnmenc.c -O pnmenc.c

echo "Step 7.4: configure, make,make install"
sleep 1
cd /usr/local/src/ffmpeg
./configure --disable-yasm
make 
make install

echo "Step 7.5: clean up"
cd /usr/local/src
rm -rf ffmpeg-snapshot.tar.bz2

echo "Done installing dependencies"
