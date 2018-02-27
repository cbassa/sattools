#!/usr/bin/env python
from astropy.coordinates import SkyCoord,EarthLocation,AltAz,get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize,interpolate



loc=EarthLocation(lat=52.8344*u.deg,lon=6.3785*u.deg,height=10*u.m)
t0=Time.now()
refalt=-6.0

# Compute solar position every 10 minutes
t=Time(t0.mjd+np.arange(144)/144.0,format='mjd',scale='utc')
psun=get_sun(t)
hor=psun.transform_to(AltAz(obstime=t,location=loc))

# Interpolating function
alt_diff=hor.alt.deg-refalt
f=interpolate.interp1d(t.mjd,alt_diff)

# Find sunrise/sunset
sign=alt_diff[1:]*alt_diff[:-1]
idx=np.where(sign<0.0)[0]
for i in idx:
    # Set
    if alt_diff[i]>alt_diff[i+1]:
        tset=Time(optimize.bisect(f,t[i].mjd,t[i+1].mjd),format='mjd',scale='utc')
    # Rise
    else:
        trise=Time(optimize.bisect(f,t[i].mjd,t[i+1].mjd),format='mjd',scale='utc')

