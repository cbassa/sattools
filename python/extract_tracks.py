#!/usr/bin/env python

import sys
from stio import fourframe,satid,observation
import numpy as np
import ppgplot
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.coordinates import SkyCoord

def pos2rel(ff,x,y):
    # Setup wcs
    w=wcs.WCS(naxis=2)
    w.wcs.crpix=ff.crpix
    w.wcs.crval=ff.crval
    w.wcs.cd=ff.cd
    w.wcs.ctype=ff.ctype
    w.wcs.set_pv([(2,1,45.0)])

    world=w.wcs_pix2world(np.array([[x,y]]),1)

    return world[0,0],world[0,1]

def dec2sex(angle):
    if angle<0.0:
        sign="-"
    else:
        sign="+"
    angle=60.0*np.abs(angle)
    

# Settings

# Minimum predicted velocity (pixels/s)
drdtmin=10.0

# Track selection region around prediction (pixels)
trkrmin=10.0

# Track selection sigma
trksig=5.0

# Minimum track points
ntrkmin=10

fname=sys.argv[1]

# Read four frame
ff=fourframe(fname)

# Read satelite IDs
try:
    f=open(fname+".id")
except OSError:
    print("Cannot open",fname+".id")
else:
    lines=f.readlines()
    f.close()

# ppgplot arrays
tr=np.array([-0.5,1.0,0.0,-0.5,0.0,1.0])
heat_l=np.array([0.0,0.2,0.4,0.6,1.0])
heat_r=np.array([0.0,0.5,1.0,1.0,1.0])
heat_g=np.array([0.0,0.0,0.5,1.0,1.0])
heat_b=np.array([0.0,0.0,0.0,0.3,1.0])

# Loop over identifications
for line in lines:
    # Decode
    id=satid(line)

    # Skip slow moving objects
    drdt=np.sqrt(id.dxdt**2+id.dydt**2)
    if drdt<drdtmin:
        continue

    # Extract significant pixels
    x,y,t,sig=ff.significant(trksig,id.x0,id.y0,id.dxdt,id.dydt,trkrmin)

    # Fit tracks
    if len(x)>ntrkmin:
        obs=observation(ff,x,y,t,sig)

        obs.ra,obs.de=pos2rel(ff,obs.x0,obs.y0)
        
        # Plot
        ppgplot.pgopen("/xs")
        ppgplot.pgpap(0.0,1.0)
        ppgplot.pgsvp(0.1,0.95,0.1,0.8)

        ppgplot.pgsch(0.8)
        ppgplot.pgmtxt("T",6.0,0.0,0.0,"UT Date: %.23s  COSPAR ID: %04d"%(ff.nfd,ff.site_id))
        ppgplot.pgmtxt("T",4.8,0.0,0.0,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')"%(ff.crval[0],3600.0*ff.crres[0],ff.crval[1],3600.0*ff.crres[1]))
        ppgplot.pgmtxt("T",3.6,0.0,0.0,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d"%(ff.wx,ff.wy,3600.0*ff.sx,3600.0*ff.sy))
        ppgplot.pgmtxt("T",2.4,0.0,0.0,"Stat: %5.1f+-%.1f (%.1f-%.1f)"%(np.mean(ff.zmax),np.std(ff.zmax),ff.vmin,ff.vmax))
        
        ppgplot.pgsch(1.0)
        ppgplot.pgwnad(0.0,ff.nx,0.0,ff.ny)
        ppgplot.pglab("x (pix)","y (pix)"," ")
        ppgplot.pgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5)
        
        ppgplot.pgimag(ff.zmax,ff.nx,ff.ny,0,ff.nx-1,0,ff.ny-1,ff.vmax,ff.vmin,tr)
        ppgplot.pgbox("BCTSNI",0.,0,"BCTSNI",0.,0)

        ppgplot.pgsci(3)
#        ppgplot.pgpt(x,y,4)

        ppgplot.pgsci(4)
        ppgplot.pgpt(np.array([obs.x0]),np.array([obs.y0]),4)
        ppgplot.pgmove(obs.xmin,obs.ymin)
        ppgplot.pgdraw(obs.xmax,obs.ymax)
        
        ppgplot.pgend()

        plt.figure()
        plt.plot(t,x,'.')
        plt.plot(t,y,'.')
        plt.plot([obs.tmin,obs.tmax],[obs.xmin,obs.xmax])
        plt.plot([obs.tmin,obs.tmax],[obs.ymin,obs.ymax])
        plt.show()

        # Track and stack
    else:
        ztrk=ff.track(id.dxdt,id.dydt,0.0)
    
#        plt.figure(figsize=(15,10))
#        plt.imshow(ztrk,origin='lower')
#        plt.scatter(id.x0,id.y0,s=150,edgecolors="y",facecolors="none")
#        plt.show()

