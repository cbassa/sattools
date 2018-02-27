#!/usr/bin/env python

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
from scipy import optimize

def residual(a,x,y,z):
    return z-(a[0]+a[1]*x+a[2]*y)

def offsets(ra0,de0,ra,de):
    w=wcs.WCS(naxis=2)
    w.wcs.crpix=np.array([0.0,0.0])
    w.wcs.crval=np.array([ra0,de0])
    w.wcs.cd=np.array([[1.0/3600.0,0.0],[0.0,1.0/3600.0]])
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    w.wcs.set_pv([(2,1,45.0)])

    world=np.stack((ra,de),axis=-1)
    pix=w.wcs_world2pix(world,1)
    return pix[:,0],pix[:,1]
    

def read_pixel_catalog(fname):
    d=np.loadtxt(fname)
    x,y,mag=d[:,0],d[:,1],d[:,2]
    
    return x,y,mag

def read_tycho2_catalog(maxmag=9.0):
    hdu=fits.open("/home/bassa/code/python/satellite/astrometry-tycho2/build/tyc2.fits")

    ra=hdu[1].data.field('RA')
    de=hdu[1].data.field('DEC')
    mag=hdu[1].data.field('MAG_VT')

    c=mag<maxmag
    print("Selected %d stars from Tycho-2 catalog with m<%.1f"%(np.sum(c),maxmag))
    return ra[c],de[c],mag[c]

def match_catalogs(ra,de,x,y,xs,ys,rmin):
    res=[]
    for i in range(len(x)):
        dx,dy=xs-x[i],ys-y[i]
        r=np.sqrt(dx*dx+dy*dy)
        if np.min(r)<rmin:
            j=np.argmin(r)
            res.append([ra[i],de[i],xs[j],ys[j]])
            
    return np.array(res)

def fit_wcs(res,x0,y0):
    ra,de,x,y=res[:,0],res[:,1],res[:,2],res[:,3]
    ra0,de0=np.mean(ra),np.mean(de)
    dx,dy=x-x0,y-y0

    for i in range(5):
        w=wcs.WCS(naxis=2)
        w.wcs.crpix=np.array([0.0,0.0])
        w.wcs.crval=np.array([ra0,de0])
        w.wcs.cd=np.array([[1.0,0.0],[0.0,1.0]])
        w.wcs.ctype=["RA---TAN","DEC--TAN"]
        w.wcs.set_pv([(2,1,45.0)])

        world=np.stack((ra,de),axis=-1)
        pix=w.wcs_world2pix(world,1)
        rx,ry=pix[:,0],pix[:,1]

        ax,cov_q,infodict,mesg,ierr=optimize.leastsq(residual,[0.0,0.0,0.0],args=(dx,dy,rx),full_output=1)
        ay,cov_q,infodict,mesg,ierr=optimize.leastsq(residual,[0.0,0.0,0.0],args=(dx,dy,ry),full_output=1)

        ra0,de0=w.wcs_pix2world([[ax[0],ay[0]]],1)[0]

    w=wcs.WCS(naxis=2)
    w.wcs.crpix=np.array([x0,y0])
    w.wcs.crval=np.array([ra0,de0])
    w.wcs.cd=np.array([[ax[1],ax[2]],[ay[1],ay[2]]])
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    w.wcs.set_pv([(2,1,45.0)])

    whdr={"CRPIX1":x0,"CRPIX2":y0,"CRVAL1":ra0,"CRVAL2":de0,"CD1_1":ax[1],"CD1_2":ax[2],"CD2_1":ay[1],"CD2_2":ay[2],"CTYPE1":"RA---TAN","CTYPE2":"DEC--TAN","CUNIT1":"DEG","CUNIT2":"DEG"}

    return whdr

# FITS file
fname="2018-02-12T17:39:55.270.fits"

# Read fits
hdu=fits.open(fname)

# Fix to force WCS with two axes only
hdu[0].header["NAXIS"]=2
w=wcs.WCS(hdu[0].header)

# Read data
data=hdu[0].data
zavg,zstd,zmax,znum=data
ny,nx=zavg.shape

# Get extrema
vmin=np.mean(zavg)-2.0*np.std(zavg)
vmax=np.mean(zavg)+3.0*np.std(zavg)

# Read pixel catalog
xs,ys,ms=read_pixel_catalog(fname+".cat")

# Read Tycho2 catalog
ra,de,mag=read_tycho2_catalog(9.0)

# Convert to pixel coordinates
xc,yc=w.all_world2pix(ra,de,1)
c=(xc>0) & (xc<nx) & (yc>0) & (yc<ny)
xt,yt=xc[c],yc[c]
rat,det=ra[c],de[c]
            
# Match catalogs
res=match_catalogs(rat,det,xt,yt,xs,ys,10.0)

whdr=fit_wcs(res,nx//2,ny//2)

    


#hdu=fits.PrimaryHDU(header=w.to_header(),data=data)
#hdr=fits.Header()
hdr=hdu[0].header
for k,v in whdr.items():
    hdr[k]=v
hdu=fits.PrimaryHDU(header=hdr,data=data)
hdu.writeto("test.fits",overwrite=True)

# Convert catalog stars
plt.figure(figsize=(10,8))
plt.imshow(zavg,origin="lower",aspect=1.0,interpolation="None",vmin=vmin,vmax=vmax)
plt.scatter(xs,ys,s=80,edgecolors="r",facecolors="none")
plt.scatter(xt,yt,s=150,edgecolors="y",facecolors="none")
#plt.scatter(xt[c],yt[c],s=50,edgecolors="k",facecolors="none")
#plt.xlim(0,720)
#plt.ylim(0,576)
plt.show()
