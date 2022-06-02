Satellite Tracking Toolkit
=========

`sattools` is a collection of tools to facilitate photographic and video satellite tracking.

Installation
------

**Requirements**

The following software and libraries are required to compile `sattools`:
* gfortran
* gcc
* make
* git
* wcslib-dev
* libgsl-dev
* pgplot5
* libpng-dev
* libx11-dev
* libjpeg-dev
* libexif-dev
* dos2unix
* sextractor

On Debian/Ubuntu, these can be installed with `sudo apt install git make dos2unix sextractor wcslib-dev pgplot5 libgsl-dev gfortran libpng-dev libx11-dev libjpeg-dev libexif-dev`.

**Build and install `qfits`**

```
wget -c ftp://ftp.eso.org/pub/qfits/qfits-5.2.0.tar.gz
gunzip -c qfits-5.2.0.tar.gz | tar xvf -
cd qfits-5.2.0
chmod +w src/xmemory.c
sed -i -e "s/swapfd = open(fname, O_RDWR | O_CREAT);/swapfd = open(fname, O_RDWR | O_CREAT, 0644);/g" src/xmemory.c
./configure
make
sudo make install
```
The `qfits` code can also be obtained from https://spacecruft.org/spacecruft/qfits

**Clone and install `sattools`**

```
git clone https://github.com/cbassa/sattools.git
cd sattools
make
```

**Setup environment variables**

You will need to set the following environment variables to run **sattools**. Add these to a login file such as `.bashrc`.
* `ST_COSPAR` COSPAR observing site number,use 9990 if you do not have one. Add your site information to `$ST_DATADIR/data/sites.txt`.
* `ST_DATADIR` path to sattools directory (e.g. `$HOME/software/sattools`)
* `ST_TLEDIR` path to TLE directory (e.g. `$HOME/tle`)
* `ST_OBSDIR` path to observations directory (e.g. `$HOME/satobs`)
* `ST_LOGIN` optional space-track.org login info for TLE download (of the form `ST_LOGIN="identity=username&password=password"`)

**Final configuration steps**

Add the `sattools` directory with executables to your `PATH` variable. On Ubuntu this is `.profile` and should be of the form `PATH=<path_to_sattools>:$PATH`.

Then reload the `.profile` file with `source $HOME/.profile`.

Try to download TLEs using `tleupdate`, which should now exist in your `PATH` variable.

**Miscellaneous setup notes**

* If you have multiple capture devices you will need to add a `/etc/udev/rules.d/99-server.rules` file to add symlinks and use them to address a particular camera. `sattools` will automatically select the camera that is scheduled for each observation.  You may use a command such as `udevadm info -a -n /dev/video0` to get your capture device attributes and use that to create the rules file. A sample rules file is available as guide in `data/`. Note that symlinks to the rules file do not work, the rules file must be modified to suit your needs and copied to `/etc/udev/rules.d/`
* You should install NTP support on the system and configure time/date to automatically synchronize to time servers.

Tools
-----

* `skymap`: Vizualize satellites on a map of the sky. 
  Example usage:
  - Compute and plot the path of a satellite (NORAD ID 43145), using a TLE from an input file (`$ST_TLEDIR/classfd.tle`) for a site (COSPAR ID 4171) at a date/time (2019-12-01T17:42:00) while looking at a particular RA/Dec (`01h00m00s`, `+30d00m00s`).
  ```
  skymap -s 4171 -t 2019-12-01T17:42:00 -c $ST_TLEDIR/classfd.tle -i 43145 -R 01:00:00 -D 30:00:00
  ```

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

* `satmap`: Visualize satellite tracks on a map of Earth
  Example usage:
  - Show the ground track of all satellites in `ST_TLEDIR/classfd.tle` (for roughly 1.3 orbits from now on)
    ```
    satmap
    ```
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
