#!/usr/bin/env python3

import argparse
from datetime import datetime, timedelta

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs84

from astropy.coordinates import Longitude, Angle
from astropy import units as u
from astropy.time import Time



def gmst(t):
    return Time(t).sidereal_time('mean', 'greenwich')

def fractional_days(t):
    d = t - datetime(t.year, 1, 1)
    return t.timetuple().tm_yday + (d.seconds + d.microseconds/1e6)/(60*60*24)

def print_tle(tle, new_epoch, new_nodeo):
    def checksum(proto_tle):
        s = 0
        for c in proto_tle:
            if c.isdigit():
                s += int(c)
            if c == '-':
                s += 1
        return s%10

    tle0,old1,old2 = tle
    tle1_proto = '{} {:2d}{:.8f} {}'.format(old1[:17],
                              abs(new_epoch.year)%100,
                              fractional_days(new_epoch),
                              old1[33:-1])
    tle2_proto = '{} {:>8.4f} {}'.format(old2[:16],
                                   Angle(new_nodeo*u.radian).degree,
                                   old2[26:68])
            
    return (tle0,
            tle1_proto + str(checksum(tle1_proto)),
            tle2_proto + str(checksum(tle2_proto)))

def launch_tle(tle, launch_date, new_launch_date):
    sat = twoline2rv(tle[1], tle[2], whichconst=wgs84)

    # New Epoch
    new_epoch = sat.epoch - launch_date + new_launch_date

    # New Right ascension of ascending node in radians
    new_nodeo = (gmst(new_launch_date) - gmst(launch_date) + Longitude(sat.nodeo * u.radian)).rad
    
    return print_tle(tle, new_epoch, new_nodeo)


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Adjust a TLE to a new launch time")
    parser.add_argument('TLEFILE', help="reference tle file", type=str)
    parser.add_argument("OLD",
                        help="Reference (old) launch time (YYYY-MM-DDTHH:MM:SS)",
                        type=str)
    parser.add_argument("NEW",
                        help="New launch time (YYYY-MM-DDTHH:MM:SS) [default: now]",
                        default=datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
                        type=str)

    args = parser.parse_args()

    with open(args.TLEFILE, 'r') as f:
        tle = f.readlines()
    tle = [l.strip() for l in tle]
 
    launch_date_ref = datetime.strptime(args.OLD, '%Y-%m-%dT%H:%M:%S')
    new_launch_date = datetime.strptime(args.NEW, '%Y-%m-%dT%H:%M:%S')
    
    new_tle = launch_tle(tle, launch_date_ref, new_launch_date)
    
    print('New TLE for a launch delay of {}h :\n'.format(new_launch_date - launch_date_ref))
    print(new_tle[0].strip()+'-delayed')
    print(new_tle[1])
    print(new_tle[2])
