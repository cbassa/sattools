Satellite Tracking Toolkit
=========

Sattools is a collection of tools to facilitate Photographic and Video satellite tracking.

Install
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

* Helper scripts install_dependencies.sh and install_sattools.sh are available at scripts directory.
  You can try run these scripts to install or use them as install guide.
  Note that install_dependencies.sh needs to be run with admin privileges (sudo ./install_dependencies.sh).

* If you re-run install_sattools.sh you should previously rmdir sattools directory or otherwise souces
  will not be fetched even if they are not present at that dir

Run notes
---------
* You will need to set the following environment variables to run **sattools**.
  * `ST_COSPAR` COSPAR number
  * `ST_DATADIR` path to sattools directory
  * `ST_TLEDIR` path to TLE directory
  * `ST_OBSDIR` path to observations directory
  * `ST_LOGIN` space-track.org login info (of the form `ST_LOGIN="identity=username&password=password"`)
  These variables are set with default values after running install_sattools.sh (except `ST_LOGIN`).
* If you have multiple capture devices you will need to add a /etc/udev/rules.d/99-server.rules file to add symlinks and use them to
  address a particular camera. Sattools will automatically select the camera that is scheduled for each observation.
  You may use a command such as 'udevadm info -a -n /dev/video0' to get your capture device attributes and
  use that to create the rules file.
  A sample rules file is available as guide in data/
  Note that symlinks to the rules file do not work, the rules file must be modified to suit your needs
  and copied to /etc/udev/rules.d/
* You should install NTP support on the system and configure time/date to automatically
  sinchronize to time servers.

Tools
-----

* `tleinfo`: Display information about a set of TLEs.
  Example usage:
  - List values (SATNO, YEAR DOY, INCL, ASCN, ARGP, MA, ECC, MM) of the TLEs in the file `bulk.tle`: `tleinfo -H -1`
  - List orbital parameters (SEMI, PERIGEE, APOGEE, PERIOD, ECC) of the TLEs in the file `bulk.tle`: `tleinfo -H -1 -a`
  - Show human-readable info of the TLE for object 74001: `tleinfo -i 74001 -f`
  
* `faketle`: Calculate a TLE based on given orbit/launch parameters.
  Example usage:
  - Assuming a standard GTO launch from Cape Canaveral (latitude 28.5°N), GTO insertion burn (10° E, 0° N) at south-bound equator crossing +1655 seconds after launch,
    launch at 2019-02-22T01:45:00, perigee/apogee heights 258/59998 km, the pre-launch TLE can be generated with
    ```
    faketle -t 2019-02-22T01:45:00 -d 1655 -q 258 -Q 59998 -I 28.5 -m 0 -w 180 -n 10
    ```

* `pass` : Calculate the next passes above a given COSPAR site.
  Example usage:
  - List all ISS passes (NORAD ID 25544) at Dwingeloo Telescope (COSPAR
    site ID 9998) for the next 10 hours (3600 seconds) based on the TLE
    in catalog.tle (drop `-a` to only show the visible passes)
    ```
    pass -i 25544 -s 9998 -l 36000 -c $ST_TLEDIR/catalog.tle -a
    ```
  - Append `-P` in order to invoke `skymap` to plot a skymap and the
    sky track for each predicted pass

* `satorbit`: Show a 3D representation of the earth, the current position,
  footprint and orbit of the selected object
  Example usage:
  - Show ISS position, footprint and orbit
    ```
    satorbit -i 25544 -c $ST_TLEDIR/catalog.tle
    ```

* `launchtle`: This tool takes an input TLE and launchtime and then corrects
  the epoch of the TLE for the new launch time, and adjusts the right ascension
  of the ascending node for how much it would have changed between the old
  and new launch time.
  Example usage:
  - Use the given pre-launch TLE of object 70002 (stored in prelaunch.tle) associated with
    the original launch time of 2018-11-11T03:00:00Z to correct for the new launch time
    at 2018-11-11T05:00:00Z and store it with a new identifier (70003 here):
    ```
    launchtle -c prelaunch.tle -i 70002 -t 2018-11-11T03:00:00 -T 2018-11-11T05:00:00 -I 70003
    ```

* `tleupdate`: Update the local database of TLEs from various
  sources (Space-Track, inttles and classfd.zip).
  Usage:
  ```
  tleupdate
  ```
