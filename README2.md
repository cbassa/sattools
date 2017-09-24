Satellite Tracking Toolkit
=========

Sattools is a collection of tools to facilitate Photographic and Video satellite tracking.

Install notes
-------------
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

* Helper scripts install_dependencies.sh and install_sattools.sh are available at scripts directory.
  You can try run these scripts to install or use them as install guide.
  Note that install_dependencies.sh needs to be run with admin privileges (sudo ./install_dependencies.sh).

* If you re-run install_sattools.sh you should previously rmdir sattools directory or otherwise souces
  will not be fetched even if they are not present at that dir

Run notes
---------
* You will need to set the following environment variables to run sattools.
	These vars are set with default values after running install_sattolls.sh.
	`ST_COSPAR` COSPAR number 
	`ST_DATADIR` path to sattools directory 
	`ST_TLEDIR` path to TLE directory
	`ST_OBSDIR` path to observations directory
* If you have multiple capture devices you will need to add a /etc/udev/rules.d/99-server.rules file to add symlinks and use them to
  address a particular camera. Sattools will automatically select the camera that is scheduled for each observation.
  You may use a command such as 'udevadm info -a -n /dev/video0' to get your capture device attributes and
  use that to create the rules file.
  A sample rules file is available as guide in data/
  Note that symlinks to the rules file do not work, the rules file must be modified to suit your needs
  and copied to /etc/udev/rules.d/
* You should install NTP support on the system and configure time/date to automatically
  sinchronize to time servers.
* Modify stget.sh for your space-track.org login and password (--post-data='identity=login&password=password')

