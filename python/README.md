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
