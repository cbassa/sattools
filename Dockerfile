# Best Practices for Writing Docker Files
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices/

# Note that this currently is a larger-than-necessary image due to the the 
# intermediate RUN commands to compile the different pages
# Stringing-it-all-together would reduce the size, but would make build
# debugging more challenging

FROM ubuntu:18.04

# To avoid questions from tzdata
ARG DEBIAN_FRONTEND=noninteractive

# Modify the following with your preferred COSPAR station ID
ENV ST_COSPAR 9999

# Update these directories per your preferences
ENV HOME       /root
ENV SAT_DIR    $HOME/satellite
ENV ST_DATADIR $SAT_DIR/sattools
ENV ST_TLEDIR  $SAT_DIR/tle
ENV ST_OBSDIR  $SAT_DIR/satobs

# Note that its convenient to mount your host filesystem to the ST_TLEDIR and ST_OBSDIR with
# -V YOUR_EXTERNAL_TLE_DIR:/root/satellite/tle
# -V YOUR_EXTERNAL_SATOBS_DIR:/root/satellite/satobs

RUN echo "Step 1.0: Install dependencies from APT" \
& apt-get update && apt-get install -y --no-install-recommends \
    git \
    make \
    dos2unix \
    sextractor \
    wcslib-dev \
    pgplot5 \
    libgsl-dev \
    gfortran \
    g++ \
    libpng-dev \
    libx11-dev \
    libjpeg-dev \
    libexif-dev \
    unzip \
    vim \
    nano \
    wget \
# For hough3dlines
    libeigen3-dev \  
# For STVID
    python3 \
    python3-dev \
    astrometry.net \
    # For scipy
    libblas3 \
    liblapack3 \
    liblapack-dev \
    libblas-dev \
# For STRF
    libfftw3-dev \
# For pip 
    ca-certificates \
&& rm -rf /var/lib/apt/lists/*

# For STVID package requirements.txt
WORKDIR /tmp
RUN ["/bin/bash", "-c", "set -o pipefail && wget -O - https://bootstrap.pypa.io/get-pip.py | python3 -"]

RUN echo "Step 2.0: Install qfits for SATTOOLS"
WORKDIR /usr/local/src
RUN ["/bin/bash", "-c", "set -o pipefail &&  wget -O - ftp://ftp.eso.org/pub/qfits/qfits-5.2.0.tar.gz | tar xzvf - \
&& cd qfits-5.2.0 \
&& chmod +w src/xmemory.c \
&& sed -i -e 's/swapfd = open(fname, O_RDWR | O_CREAT);/swapfd = open(fname, O_RDWR | O_CREAT, 0644);/g' src/xmemory.c \
&& ./configure \
&& make \
&& make install \
&& echo 'Step 2.1: clean up' \
&& rm -rf /usr/local/src/qfits-5.2.0"]

WORKDIR ${SAT_DIR}
RUN echo "Step 3.0: Install hough3dlines for STVID" \
&& git clone https://gitlab.com/pierros/hough3d-code.git \
&& cd hough3d-code \
&& make \
&& rm *.o \
&& echo "Done installing dependencies" \
&& mkdir -p $SAT_DIR

##
## THIS SECTION BUILDS SATTOOLS SUITE
##

WORKDIR $SAT_DIR
RUN echo "Step 4.0: Build satools/strf/stvid" \
&& echo "Step 4.1: Build sattools: Satellite Tracking Toolkit" \
&& git clone https://github.com/interplanetarychris/sattools.git \
&& cd sattools \
&& git checkout Docker \
&& make \
&& make clean

RUN echo "Step 4.2. Build stvid: Satellite tracking with video cameras" \
&& git clone https://github.com/interplanetarychris/stvid.git \
&& cd $SAT_DIR/stvid \
&& git checkout Docker \
&& pip install -r $SAT_DIR/stvid/requirements.txt \
# Needed twice, because ppgplot installs correctly on the second try
&& pip install -r $SAT_DIR/stvid/requirements.txt
#WORKDIR /usr/local/src
#RUN git clone https://github.com/haavee/ppgplot.git
# WORKDIR /usr/local/src/ppgplot
# RUN python3 setup.py install

WORKDIR $SAT_DIR
RUN echo "Step 4.3. Build strf: Radio Frequency Satellite Tracking" \
&& git clone https://github.com/cbassa/strf.git \
&& cd $SAT_DIR/strf \
&& make \
&& make clean

RUN echo "Step 4.4 Decrease container size after everything is compiled" \
&& apt-get purge -y \
    wcslib-dev \
    libgsl-dev \
    gfortran \
    libpng-dev \
    libx11-dev \
    libjpeg-dev \
    libexif-dev \
# For hough3dlines
    libeigen3-dev \  
# For STVID
    python3-dev \
    # For scipy
    liblapack-dev \
    libblas-dev \
# For STRF
    libfftw3-dev \
&& apt-get autoremove -y

WORKDIR $HOME
RUN echo "Step 5.0: set environment variables" \
&& echo "export ST_COSPAR=$ST_COSPAR"   >>$HOME/.profile \
&& echo "export ST_DATADIR=$ST_DATADIR" >>$HOME/.profile \
&& echo "export ST_TLEDIR=$ST_TLEDIR"   >>$HOME/.profile \
&& echo "export ST_OBSDIR=$ST_OBSDIR"   >>$HOME/.profile \
&& mkdir -p $ST_OBSDIR \
&& mkdir $ST_OBSDIR/control \
&& echo "Step 5.1: set path" \
&& echo "PATH=$ST_DATADIR:$ST_DATADIR/scripts:$SAT_DIR/hough3d-code:$SAT_DIR/strf:\$PATH" >>$HOME/.profile 

RUN echo "Step 5.2: Download initial TLEs" \
&& mkdir -p $ST_TLEDIR \
&& $ST_DATADIR/tleupdate \
&& echo "Final step: run" \
&& echo "source $HOME/.profile"

# Run bash by default if everything else compiles
ENTRYPOINT /bin/bash
CMD bash