## Notes


## Setup
- Dependency installation via pip:
  ```
  pip install -r requirements.txt
  ```

- In order to be able to use `ground_track.py`,
  download [earth.png](https://raw.githubusercontent.com/galactics/beyond/master/doc/source/_static/earth.png) into this directory.

## Usage

- Plot the ground track of satellite based on it's TLE around the epoch:
  ```
  ./ground_track.py ../examples/sathyabamasat.txt
  ```

- Calculate and plot satellite passes (requires a site.txt)
  ```
  ./pass.py ../examples/sathyabamasat.txt $ST_DATADIR/sites.txt -s 7300 --starttime 2019-11-06T19:30:00
  ```

- Adjust a TLE to a new launch time
  ```
  # Data from https://community.libre.space/t/electron-its-business-time-launch-this-weekend-irvine-01-amateur-payload/2819/4
  $ ./launchtle.py ../examples/irvine01.txt 2018-11-11T03:00:00 2018-11-11T04:05:00
  New TLE for a launch delay of 1:05:00h :

  IRVINE01-delayed
  1 70002U 18599A   18315.20665747  .00000000  00000-0  00000-0 0    09
  2 70002  85.1205 106.4513 0012705 292.5520 107.9249 15.20792276    05
  ```
