#!/usr/bin/env python

import sys,os,glob
from stio import fourframe,satid,observation
import numpy as np
import ppgplot
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.coordinates import SkyCoord

# Get COSPAR ID
def get_cospar(norad):
    f=open(os.getenv("ST_DATADIR")+"/data/desig.txt")
    lines=f.readlines()
    f.close()

    cospar=([line for line in lines if str(norad) in line])[0].split()[1]
    
    return "%2s %s"%(cospar[0:2],cospar[2:])

# IOD position format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
def format_position(ra,de):
    ram=60.0*ra/15.0
    rah=int(np.floor(ram/60.0))
    ram-=60.0*rah
    
    des=np.sign(de)
    dem=60.0*np.abs(de)
    ded=int(np.floor(dem/60.0))
    dem-=60.0*ded

    if des==-1:
        sign="-"
    else:
        sign="+"

    return ("%02d%06.3f%c%02d%05.2f"%(rah,ram,sign,ded,dem)).replace(".","")

# Format IOD line
def format_iod_line(norad,cospar,site_id,t,ra,de):
    pstr=format_position(ra,de)
    tstr=t.replace("-","").replace("T","").replace(":","").replace(".","")

    return "%05d %-9s %04d G %s 17 25 %s 37 S"%(norad,cospar,site_id,tstr,pstr)

# Extract tracks
def extract_tracks(fname,trkrmin,drdtmin,trksig,ntrkmin):
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
        if len(t)>ntrkmin:
            obs=observation(ff,x,y,t,sig)
            cospar=get_cospar(id.norad)
            iod_line="%s"%format_iod_line(id.norad,cospar,ff.site_id,obs.nfd,obs.ra,obs.de)

            print(iod_line)

            if id.catalog.find("classfd.tle")>0:
                outfname="classfd.dat"
            elif id.catalog.find("inttles.tle")>0:
                outfname="inttles.dat"
            else:
                outfname="catalog.dat"

            f=open(outfname,"a")
            f.write("%s\n"%iod_line);
            f.close()
            
            # Plot
            ppgplot.pgopen(fname.replace(".fits","")+"_%05d.png/png"%id.norad)
            ppgplot.pgpap(0.0,1.0)
            ppgplot.pgsvp(0.1,0.95,0.1,0.8)
            
            ppgplot.pgsch(0.8)
            ppgplot.pgmtxt("T",6.0,0.0,0.0,"UT Date: %.23s  COSPAR ID: %04d"%(ff.nfd,ff.site_id))
            ppgplot.pgmtxt("T",4.8,0.0,0.0,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')"%(ff.crval[0],3600.0*ff.crres[0],ff.crval[1],3600.0*ff.crres[1]))
            ppgplot.pgmtxt("T",3.6,0.0,0.0,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d"%(ff.wx,ff.wy,3600.0*ff.sx,3600.0*ff.sy))
            ppgplot.pgmtxt("T",2.4,0.0,0.0,"Stat: %5.1f+-%.1f (%.1f-%.1f)"%(np.mean(ff.zmax),np.std(ff.zmax),ff.vmin,ff.vmax))
            ppgplot.pgmtxt("T",0.3,0.0,0.0,iod_line)
        
            ppgplot.pgsch(1.0)
            ppgplot.pgwnad(0.0,ff.nx,0.0,ff.ny)
            ppgplot.pglab("x (pix)","y (pix)"," ")
            ppgplot.pgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5)
            
            ppgplot.pgimag(ff.zmax,ff.nx,ff.ny,0,ff.nx-1,0,ff.ny-1,ff.vmax,ff.vmin,tr)
            ppgplot.pgbox("BCTSNI",0.,0,"BCTSNI",0.,0)
            ppgplot.pgstbg(1)

            ppgplot.pgsci(0)
            if id.catalog.find("classfd.tle")>0:
                ppgplot.pgsci(4)
            elif id.catalog.find("inttles.tle")>0:
                ppgplot.pgsci(3)
            ppgplot.pgpt(np.array([obs.x0]),np.array([obs.y0]),4)
            ppgplot.pgmove(obs.xmin,obs.ymin)
            ppgplot.pgdraw(obs.xmax,obs.ymax)
            ppgplot.pgsch(0.65)
            ppgplot.pgtext(np.array([obs.x0]),np.array([obs.y0])," %05d"%id.norad)
            ppgplot.pgsch(1.0)
            ppgplot.pgsci(1)
            
            ppgplot.pgend()


            # Track and stack
#        else:
#            ztrk=ff.track(id.dxdt,id.dydt,0.0)
            
#            plt.figure(figsize=(15,10))
#            plt.imshow(ztrk,origin='lower')
#            plt.scatter(id.x0,id.y0,s=150,edgecolors="y",facecolors="none")
#            plt.show()


if __name__ == '__main__':
    # Minimum predicted velocity (pixels/s)
    drdtmin=10.0

    # Track selection region around prediction (pixels)
    trkrmin=10.0

    # Track selection sigma
    trksig=5.0

    # Minimum track points
    ntrkmin=10

#    extract_tracks("2018-02-26T05:26:15.801.fits",trkrmin,drdtmin,trksig,ntrkmin)
    
    files=sorted(glob.glob("2018-03-04T02*.fits"))

    for file in files:
        print(file)
        extract_tracks(file,trkrmin,drdtmin,trksig,ntrkmin)


