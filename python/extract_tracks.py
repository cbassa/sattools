#!/usr/bin/env python

import sys,os,glob
from stio import fourframe,satid,observation
import numpy as np
import ppgplot
import matplotlib.pyplot as plt
from scipy import optimize,ndimage
from astropy import wcs
from astropy.coordinates import SkyCoord

# Gaussian model
def model(a,nx,ny):
    x,y=np.meshgrid(np.arange(nx),np.arange(ny))
    dx,dy=(x-a[0])/a[2],(y-a[1])/a[2]
    arg=-0.5*(dx**2+dy**2)
    return a[3]*np.exp(arg)+a[4]

# Residual function
def residual(a,img):
    ny,nx=img.shape
    mod=model(a,nx,ny)
    return (img-mod).ravel()

# Find peak
def peakfind(img,w=1.0):
    # Find approximate location
    ny,nx=img.shape
    i=np.argmax(img)
    y0=int(i/nx)
    x0=i-y0*nx

    # Image properties
    imgavg=np.mean(img)
    imgstd=np.std(img)
    
    # Estimate
    a=np.array([x0,y0,w,img[y0,x0]-imgavg,imgavg])
    q,cov_q,infodict,mesg,ier=optimize.leastsq(residual,a,args=(img),full_output=1)

    # Extract
    xc,yc,w=q[0],q[1],q[2]

    # Significance
    sigma=(a[3]-imgavg)/imgstd

    return xc,yc,w,sigma

# Plot selection
def plot_selection(id,x0,y0,dt=2.0,w=10.0):
    dx,dy=id.x1-id.x0,id.y1-id.y0
    ang=np.arctan2(dy,dx)
    r=np.sqrt(dx**2+dy**2)
    drdt=r/(id.t1-id.t0)
    sa,ca=np.sin(ang),np.cos(ang)

    dx=np.array([-dt,-dt,dt,dt,-dt])*drdt
    dy=np.array([w,-w,-w,w,w])
    x=ca*dx-sa*dy+x0
    y=sa*dx+ca*dy+y0

    ppgplot.pgsci(7)
    ppgplot.pgline(x,y)
    
    return

# Check if point is inside selection
def inside_selection(id,x0,y0,xmid,ymid,dt=2.0,w=10.0):
    dx,dy=id.x1-id.x0,id.y1-id.y0
    ang=-np.arctan2(dy,dx)
    r=np.sqrt(dx**2+dy**2)
    drdt=r/(id.t1-id.t0)
    sa,ca=np.sin(ang),np.cos(ang)

    dx,dy=x0-xmid,y0-ymid
    rm=ca*dx-sa*dy
    wm=sa*dx+ca*dy
    dtm=rm/drdt
    
    if (abs(wm)<w) & (abs(dtm)<dt):
        return True
    else:
        return False

# Get COSPAR ID
def get_cospar(norad):
    f=open(os.getenv("ST_DATADIR")+"/data/desig.txt")
    lines=f.readlines()
    f.close()

    try:
        cospar=([line for line in lines if str(norad) in line])[0].split()[1]
    except IndexError:
        cospar="18500A"
    
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
            # Get times
            tmin=np.min(t)
            tmax=np.max(t)
            tmid=0.5*(tmax+tmin)
            mjd=ff.mjd+tmid/86400.0

            # Very simple polynomial fit; no weighting, no cleaning
            px=np.polyfit(t-tmid,x,1)
            py=np.polyfit(t-tmid,y,1)

            # Extract results
            x0,y0=px[1],py[1]
            dxdt,dydt=px[0],py[0]
            xmin=x0+dxdt*(tmin-tmid)
            ymin=y0+dydt*(tmin-tmid)
            xmax=x0+dxdt*(tmax-tmid)
            ymax=y0+dydt*(tmax-tmid)

            cospar=get_cospar(id.norad)
            obs=observation(ff,mjd,x0,y0)
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
            #ppgplot.pgopen("/xs")
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
            ppgplot.pgpt(np.array([x0]),np.array([y0]),4)
            ppgplot.pgmove(xmin,ymin)
            ppgplot.pgdraw(xmax,ymax)
            ppgplot.pgsch(0.65)
            ppgplot.pgtext(np.array([x0]),np.array([y0])," %05d"%id.norad)
            ppgplot.pgsch(1.0)
            ppgplot.pgsci(1)
            
            ppgplot.pgend()


        elif id.catalog.find("classfd.tle")>0:
            # Track and stack
            t=np.linspace(0.0,ff.texp)
            x,y=id.x0+id.dxdt*t,id.y0+id.dydt*t
            c=(x>0) & (x<ff.nx) & (y>0) & (y<ff.ny)

            # Skip if no points selected
            if np.sum(c)==0:
                continue

            # Compute track
            tmid=np.mean(t[c])
            mjd=ff.mjd+tmid/86400.0
            xmid=id.x0+id.dxdt*tmid
            ymid=id.y0+id.dydt*tmid
            ztrk=ndimage.gaussian_filter(ff.track(id.dxdt,id.dydt,tmid),1.0)
            vmin=np.mean(ztrk)-2.0*np.std(ztrk)
            vmax=np.mean(ztrk)+6.0*np.std(ztrk)

            # Select region
            xmin=int(xmid-100)
            xmax=int(xmid+100)
            ymin=int(ymid-100)
            ymax=int(ymid+100)
            if xmin<0: xmin=0
            if ymin<0: ymin=0
            if xmax>ff.nx: xmax=ff.nx-1
            if ymax>ff.ny: ymax=ff.ny-1

            # Find peak
            x0,y0,w,sigma=peakfind(ztrk[ymin:ymax,xmin:xmax])
            x0+=xmin
            y0+=ymin

            # Skip if peak is not significant
            if sigma<trksig:
                continue

            # Skip if point is outside selection area
            if inside_selection(id,xmid,ymid,x0,y0)==False:
                continue;
            
            # Format IOD line
            cospar=get_cospar(id.norad)
            obs=observation(ff,mjd,x0,y0)
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

            ppgplot.pgimag(ztrk,ff.nx,ff.ny,0,ff.nx-1,0,ff.ny-1,vmax,vmin,tr)
            ppgplot.pgbox("BCTSNI",0.,0,"BCTSNI",0.,0)
            ppgplot.pgstbg(1)

            plot_selection(id,xmid,ymid)
            
            ppgplot.pgsci(0)
            if id.catalog.find("classfd.tle")>0:
                ppgplot.pgsci(4)
            elif id.catalog.find("inttles.tle")>0:
                ppgplot.pgsci(3)
            ppgplot.pgpt(np.array([id.x0]),np.array([id.y0]),17)
            ppgplot.pgmove(id.x0,id.y0)
            ppgplot.pgdraw(id.x1,id.y1)
            ppgplot.pgpt(np.array([x0]),np.array([y0]),4)
            ppgplot.pgsch(0.65)
            ppgplot.pgtext(np.array([id.x0]),np.array([id.y0])," %05d"%id.norad)
            ppgplot.pgsch(1.0)
            ppgplot.pgsci(1)
            
            ppgplot.pgend()


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
    
    files=sorted(glob.glob("2*.fits"))

    for file in files:
        extract_tracks(file,trkrmin,drdtmin,trksig,ntrkmin)


