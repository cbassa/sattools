#!/usr/bin/env python3

"""
Script showing the position of the ISS at the time of the TLE
and the ground track for the previous and the next orbit

(original) source: https://github.com/galactics/beyond/blob/master/doc/source/_static/ground-track.py
"""

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

from beyond.io.tle import Tle
from beyond.dates import Date, timedelta


def tle_from_file(filename):
    '''
    Returns: TLE:Tle, Object name:str
    '''

    with open(filename, 'r') as f:
        lines = f.readlines()

    return Tle(''.join(lines)), lines[0].strip()


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Plot the ground track of a TLE")
    parser.add_argument('TLEFILE', help="the tle file", type=str)

    args = parser.parse_args()

    # Parsing of TLE
    tle, name = tle_from_file(args.TLEFILE)
   
    # Conversion into `Orbit` object
    orb = tle.orbit()
    
    # Tables containing the positions of the ground track
    latitudes, longitudes = [], []
    prev_lon, prev_lat = None, None
    
    period = orb.infos.period
    start = orb.date - period
    stop = 2 * period
    step = period / 100
    
    for point in orb.ephemeris(start=start, stop=stop, step=step):
    
        # Conversion to earth rotating frame
        point.frame = 'ITRF'
    
        # Conversion from cartesian to spherical coordinates (range, latitude, longitude)
        point.form = 'spherical'
    
        # Conversion from radians to degrees
        lon, lat = np.degrees(point[1:3])
    
        # Creation of multiple segments in order to not have a ground track
        # doing impossible paths
        if prev_lon is None:
            lons = []
            lats = []
            longitudes.append(lons)
            latitudes.append(lats)
        elif orb.i < np.pi /2 and (np.sign(prev_lon) == 1 and np.sign(lon) == -1):
            lons.append(lon + 360)
            lats.append(lat)
            lons = [prev_lon - 360]
            lats = [prev_lat]
            longitudes.append(lons)
            latitudes.append(lats)
        elif orb.i > np.pi/2 and (np.sign(prev_lon) == -1 and np.sign(lon) == 1):
            lons.append(lon - 360)
            lats.append(lat)
            lons = [prev_lon + 360]
            lats = [prev_lat]
            longitudes.append(lons)
            latitudes.append(lats)
    
        lons.append(lon)
        lats.append(lat)
        prev_lon = lon
        prev_lat = lat
    

    fig = plt.figure(figsize=(15.2, 8.2))
    ax = fig.add_subplot(111)

    ax.set_title(name)

    # Plot earth
    img = Path(__file__).parent / "earth.png"
    im = plt.imread(str(img))
    ax.imshow(im, extent=[-180, 180, -90, 90])

    # Plot ground track 
    for lons, lats in zip(longitudes, latitudes):
        ax.plot(lons, lats, 'r-')

    # Plot location at epoch time 
    lon, lat = np.degrees(orb.copy(frame='ITRF', form='spherical')[1:3])
    ax.plot([lon], [lat], 'bo')

    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])
    ax.grid(True, color='w', linestyle=":", alpha=0.4)
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    fig.tight_layout()
    
    if "no-display" not in sys.argv:
        plt.show()
