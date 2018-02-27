#!/usr/bin/env python
from skyfield.positionlib import ICRF
from skyfield.api import load,Topos

ts=load.timescale()
t=ts.now()
observer=Topos(latitude_degrees=52.8344,longitude_degrees=6.3785,elevation_m=10.0)
p=ICRF([0.0,0.0,0.0],observer_data=observer,t=t)
q=p.from_altaz(alt_degrees=30.0,az_degrees=30.0)

print(p,q)
