FROM ubuntu:16.04

RUN apt-get update && apt-get upgrade -y && apt-get install -y ntp eog emacs gfortran libpng-dev libx11-dev libjpeg-dev libexif-dev git dos2unix sextractor
RUN apt-get install -y build-essential vim wget

RUN cd /usr/local/src && git clone https://github.com/cbassa/sattools.git

RUN apt-get install -y pgplot5

# qfits
RUN cd /usr/local/src && wget -c ftp://ftp.eso.org/pub/qfits/qfits-5.2.0.tar.gz && tar xzf qfits-5.2.0.tar.gz
RUN cd /usr/local/src/qfits-5.2.0 && sed -i -e "s/swapfd = open(fname, O_RDWR | O_CREAT);/swapfd = open(fname, O_RDWR | O_CREAT, 0644);/g" src/xmemory.c
RUN cd /usr/local/src/qfits-5.2.0 && ./configure && make && make install

# wcslib 2.9
RUN cd /usr/local/src && wget -q -c "https://drive.google.com/uc?export=download&id=0B-15JZVdjJi4QW0zZmZUM1ZXblU" -O wcslib-2.9.tar && tar xf wcslib-2.9.tar
RUN cd /usr/local/src/wcslib-2.9/C && make clean && rm libwcs_c.a && make
RUN ln -s /usr/local/src/wcslib-2.9/C/libwcs_c.a /usr/local/lib/ && \
    ln -s /usr/local/src/wcslib-2.9/C/proj.h /usr/local/include/ && \
    ln -s /usr/local/src/wcslib-2.9/C/cel.h /usr/local/include/

# GSL 1.16
RUN cd /usr/local/src && wget -q -c ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz && tar xzf gsl-1.16.tar.gz
RUN cd /usr/local/src/gsl-1.16 && ./configure && make && make install

# FFTW3
RUN apt-get install -y fftw3-dev

# ffmpeg
RUN cd /usr/local/src && wget -q -c http://ffmpeg.org/releases/ffmpeg-snapshot.tar.bz2 && tar xjf ffmpeg-snapshot.tar.bz2
RUN cd /usr/local/src/ffmpeg/libavcodec && cat /usr/local/src/sattools/pnmenc.patch | patch -p0
RUN cd /usr/local/src/ffmpeg && ./configure --disable-yasm && make && make install

RUN cd /usr/local/src/sattools && make
RUN apt-get install -y unzip
RUN mkdir -p /usr/local/src/sattools/tle && cd /usr/local/src/sattools/tle && wget -q -c https://www.prismnet.com/~mmccants/tles/classfd.zip && unzip classfd.zip && dos2unix classfd.tle

ENV ST_COSPAR 4171
ENV ST_DATADIR /usr/local/src/sattools
ENV ST_TLEDIR /usr/local/src/sattools/tle
ENV ST_OBSDIR $HOME
