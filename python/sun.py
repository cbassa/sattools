#!/usr/bin/env python
from astropy.coordinates import SkyCoord,EarthLocation,AltAz,get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np
from scipy import optimize,interpolate
import time


def get_sunriseset(t,lat,lon,height,refalt):
    loc=EarthLocation(lat=lat*u.deg,lon=lon*u.deg,height=height*u.m)
    t0=Time(t,format='unix')
    print(t0.isot)
    # Compute solar position every 10 minutes
    t=Time(t0.mjd+np.arange(160)/160.0,format='mjd',scale='utc')
    psun=get_sun(t)
    hor=psun.transform_to(AltAz(obstime=t,location=loc))
    
    # Interpolating function
    alt_diff=hor.alt.deg-refalt
    f=interpolate.interp1d(t.mjd,alt_diff)
    
    # Find sunrise/sunset
    sign=alt_diff[1:]*alt_diff[:-1]
    idx=np.where(sign<0.0)[0]
    print(idx)
    for i in idx:
        # Set
        if alt_diff[i]>alt_diff[i+1]:
            tset=Time(optimize.bisect(f,t[i].mjd,t[i+1].mjd),format='mjd',scale='utc')
            # Rise
        else:
            trise=Time(optimize.bisect(f,t[i].mjd,t[i+1].mjd),format='mjd',scale='utc')

    return trise.unix,tset.unix

if __name__ == '__main__':
    lat,lon,height=52.8344,6.3785,10

    t=float(time.time())
    t=1520573400.0
    print(t,type(t))
    trise,tset=get_sunriseset(t,lat,lon,height,-6.0)

    print(t,trise,tset)
