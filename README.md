Satellite Tracking Toolkit
=========

Sattools is a collection of tools to facilitate Photographic and Video satellite tracking.

Build
------
* Clone locally the code repository
* Install common dependencies
  * gfortran
  * gcc
  * libpng-dev
  * libx11-dev
  * libjpeg-dev 
  * libexif-dev
* Build & install required libraries
  * qfits-5.2.0: ftp://ftp.eso.org/pub/qfits/qfits-5.2.0.tar.gz
  * pgplot-5.2.2: http://www.astro.caltech.edu/~tjp/pgplot/
  * gsl-1.15: ftp://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz 
  * wcslib-2.9: http://www.epta.eu.org/~bassa/wcslib-2.9.tar
* Run `make` on the sattools folder

Run
-----
You will need to set the following environment variables to run sattools:
* `ST_COSPAR` COSPAR number 
* `ST_DATADIR` path to sattools directory 
* `ST_TLEDIR` path to TLE directory
