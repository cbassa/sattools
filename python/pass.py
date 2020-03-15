#!/usr/bin/env python

"""
Script to calculate the next passes of a given satellite (from TLE) at a given location
"""
import sys
import argparse
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u

from beyond.dates import Date, timedelta
from beyond.env.jpl import get_orbit
from beyond.io.tle import Tle
from beyond.frames import create_station
from beyond.config import config

from utils import tle_from_file, read_site_file


def polar_plot(pass_data):
    plt.figure()
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    plt.plot(np.radians(pass_data['azims']), pass_data['elevs'], '-')

    pass_data['meta'] = {}
    for i, event in enumerate(pass_data['event']):
        if event:
            if event.info == 'AOS':
                style = 'ro'
                pass_data['meta'].update({'AOS': i})
            elif event.info == 'MAX':
                style = 'go'
                pass_data['meta'].update({'MAX': i})
            elif event.info == 'LOS':
                style = 'bo'
                pass_data['meta'].update({'LOS': i})
            else:
                style = 'bo'
            plt.plot(np.radians(pass_data['azims'][i]), pass_data['elevs'][i], style)

    ax.set_yticks(range(0, 90, 20))
    ax.set_yticklabels(map(str, range(90, 0, -20)))
    ax.set_rmax(90)

    try:
        aos_str = pass_data['date'][pass_data['meta']['AOS']].datetime.strftime('%H:%M:%S')
        max_str = pass_data['date'][pass_data['meta']['MAX']].datetime.strftime('%H:%M:%S')
        los_str = pass_data['date'][pass_data['meta']['LOS']].datetime.strftime('%H:%M:%S')
        textstr = 'AOS {}\nMAX {}\n LOS {}'.format(aos_str, max_str, los_str)

        ax.text(0.05,
                0.95,
                textstr,
                transform=ax.transAxes, 
                fontsize=14,
                verticalalignment='top',
                bbox={'boxstyle': 'round',
                      'facecolor': 'white',
                      'alpha': 0.5})
    except KeyError:
        # TODO: Fix partial passes
        pass

    plt.show()


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Compute the next passes of a satellite at a given location")
    parser.add_argument('TLEFILE', help="the tle file", type=str)
    parser.add_argument('SITEFILE', help="the site catalog file", type=str)
    parser.add_argument("-s",
                        "--site",
                        help="COSPAR_ID (must be present in sites.txt",
                        type=int)
    parser.add_argument("-t",
                        "--starttime",
                        help="Start time (YYYY-MM-DDTHH:MM:SS) [default: now]",
                        default=datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"))
    parser.add_argument("-d",
                        "--duration",
                        help="Duration to search for passes [hours; default: 1.0]",
                        type=float,
                        default=1.0)
    args = parser.parse_args()

    sites = read_site_file(args.SITEFILE)

    if not args.site in sites.keys():
        print("ERROR: Site {} not found in {}.".format(args.site, args.SITEFILE))
        sys.exit(-1)

    site = sites[args.site]

    tle, name = tle_from_file(args.TLEFILE)
    orbit = tle.orbit()

    # Station definition
    station = create_station(site['symbol'], (site['latitude'].to(u.deg).value,
                                              site['longitude'].to(u.deg).value,
                                              site['altitude'].to(u.m).value))
    starttime = datetime.strptime(args.starttime, '%Y-%m-%dT%H:%M:%S')

    pass_data = {'azims': [],
                 'elevs': [],
                 'date': [],
                 'r': [],
                 'r_dot': [],
                 'event': []}

    print("    Time      Azim    Elev    Distance   Radial Velocity")
    print("=========================================================")
    
    for orb in station.visibility(orbit,
                                  start=Date(starttime),
                                  stop=timedelta(hours=args.duration),
                                  step=timedelta(seconds=30),
                                  events=True):
        elev = np.degrees(orb.phi)
        # Radians are counterclockwise and azimuth is clockwise
        azim = np.degrees(-orb.theta) % 360
    
        # Archive for plotting
        pass_data['azims'].append(azim)
        # Matplotlib actually force 0 to be at the center of the polar plot,
        # so we trick it by inverting the values
        pass_data['elevs'].append(90 - elev)
        pass_data['date'].append(orb.date)
        pass_data['r'].append(orb.r)
        pass_data['r_dot'].append(orb.r_dot)
        pass_data['event'].append(orb.event)

        r = orb.r / 1000.
        print("{event:7} {orb.date:%H:%M:%S} {azim:7.2f} {elev:7.2f} {r:10.2f} {orb.r_dot:10.2f}".format(
            orb=orb, r=r, azim=azim, elev=elev, event=orb.event.info if orb.event is not None else ""
        ))
    
        if orb.event and orb.event.info.startswith("LOS"):
            # End of the pass, plot the result and
            # clear the arrays for the next pass
            print()
            polar_plot(pass_data)
            pass_data = {'azims': [],
                         'elevs': [],
                         'date': [],
                         'r': [],
                         'r_dot': [],
                         'event': []}

    # Plot the last (incomplete) pass 
    if pass_data['azims']:
        polar_plot(pass_data)
