#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy import wcs

class observation:
    """Satellite observation"""

    def __init__(self,ff,mjd,x0,y0):
        """Define an observation"""

        # Store
        self.mjd=mjd
        self.x0=x0
        self.y0=y0
        
        # Get times
        self.nfd=Time(self.mjd,format='mjd',scale='utc').isot
        
        # Correct for rotation
        tobs=Time(ff.mjd+0.5*ff.texp/86400.0,format='mjd',scale='utc')
        tobs.delta_ut1_utc=0
        hobs=tobs.sidereal_time("mean",longitude=0.0).degree
        tmid=Time(self.mjd,format='mjd',scale='utc')
        tmid.delta_ut1_utc=0
        hmid=tmid.sidereal_time("mean",longitude=0.0).degree
        
        # Compute ra/dec
        world=ff.w.wcs_pix2world(np.array([[self.x0,self.y0]]),1)
        self.ra=world[0,0]+hobs-hmid
        self.de=world[0,1]

class satid:
    """Satellite identifications"""

    def __init__(self,line):
        s=line.split()
        self.nfd=s[0]
        self.x0=float(s[1])
        self.y0=float(s[2])
        self.t0=0.0
        self.x1=float(s[3])
        self.y1=float(s[4])
        self.t1=float(s[5])
        self.norad=int(s[6])
        self.catalog=s[7]
        self.state=s[8]
        self.dxdt=(self.x1-self.x0)/(self.t1-self.t0)
        self.dydt=(self.y1-self.y0)/(self.t1-self.t0)

    def __repr__(self):
        return "%s %f %f %f -> %f %f %f %d %s %s"%(self.nfd,self.x0,self.y0,self.t0,self.x1,self.y1,self.t1,self.norad,self.catalog,self.state)
        
class fourframe:
    """Four frame class"""

    def __init__(self,fname=None):
        if fname==None:
            # Initialize empty fourframe
            self.nx=0
            self.ny=0
            self.nz=0
            self.mjd=-1
            self.nfd=None
            self.zavg=None
            self.zstd=None
            self.zmax=None
            self.znum=None
            self.dt=None
            self.site_id=0
            self.observer=None
            self.texp=0.0
            self.fname=None
            self.crpix=np.array([0.0,0.0])
            self.crval=np.array([0.0,0.0])
            self.cd=np.array([[1.0,0.0],[0.0,1.0]])
            self.ctype=["RA---TAN","DEC--TAN"]
            self.cunit=np.array(["deg","deg"])
            self.crres=np.array([0.0,0.0])
        else:
            # Read FITS file
            hdu=fits.open(fname)

            # Read image planes
            self.zavg,self.zstd,self.zmax,self.znum=hdu[0].data

            # Frame properties
            self.ny,self.nx=self.zavg.shape
            self.nz=hdu[0].header['NFRAMES']
    
            # Read frame time oselfsets
            self.dt=np.array([hdu[0].header['DT%04d'%i] for i in range(self.nz)])
    
            # Read header
            self.mjd=hdu[0].header['MJD-OBS']
            self.nfd=hdu[0].header['DATE-OBS']
            self.site_id=hdu[0].header['COSPAR']
            self.observer=hdu[0].header['OBSERVER']
            self.texp=hdu[0].header['EXPTIME']
            self.fname=fname

            # Astrometry keywords
            self.crpix=np.array([hdu[0].header['CRPIX1'],hdu[0].header['CRPIX2']])
            self.crval=np.array([hdu[0].header['CRVAL1'],hdu[0].header['CRVAL2']])
            self.cd=np.array([[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],
                              [hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]])
            self.ctype=[hdu[0].header['CTYPE1'],hdu[0].header['CTYPE2']]
            self.cunit=[hdu[0].header['CUNIT1'],hdu[0].header['CUNIT2']]
            self.crres=np.array([hdu[0].header['CRRES1'],hdu[0].header['CRRES2']])

        # Compute image properties
        self.sx=np.sqrt(self.cd[0,0]**2+self.cd[1,0]**2)
        self.sy=np.sqrt(self.cd[0,1]**2+self.cd[1,1]**2)
        self.wx=self.nx*self.sx
        self.wy=self.ny*self.sy
        self.vmin=np.mean(self.zmax)-2.0*np.std(self.zmax)
        self.vmax=np.mean(self.zmax)+6.0*np.std(self.zmax)

        # Setup WCS
        self.w=wcs.WCS(naxis=2)
        self.w.wcs.crpix=self.crpix
        self.w.wcs.crval=self.crval
        self.w.wcs.cd=self.cd
        self.w.wcs.ctype=self.ctype
        self.w.wcs.set_pv([(2,1,45.0)])
        

    def significant(self,sigma,x0,y0,dxdt,dydt,rmin=10.0):
        """Extract significant points"""

        # Generate sigma frame
        zsig=(self.zmax-self.zavg)/(self.zstd+1e-9)

        # Select
        c=(zsig>sigma)

        # Positions
        xm,ym=np.meshgrid(np.arange(self.nx),np.arange(self.ny))
        x,y=np.ravel(xm[c]),np.ravel(ym[c])
        inum=np.ravel(self.znum[c]).astype('int')
        sig=np.ravel(zsig[c])
        t=np.array([self.dt[i] for i in inum])

        # Predicted positions
        xr=x0+dxdt*t
        yr=y0+dydt*t
        r=np.sqrt((x-xr)**2+(y-yr)**2)
        c=r<rmin

        return x[c],y[c],t[c],sig[c]

    def track(self,dxdt,dydt,tref):
        """Track and stack"""
        # Empty frame
        ztrk=np.zeros_like(self.zavg)
        
        # Loop over frames
        for i in range(self.nz):
            dx=int(np.round(dxdt*(self.dt[i]-tref)))
            dy=int(np.round(dydt*(self.dt[i]-tref)))

            # Skip if shift larger than image
            if np.abs(dx)>=self.nx:
                continue
            if np.abs(dy)>=self.ny:
                continue

            # Extract range
            if dx>=0:
                i1min,i1max=dx,self.nx-1
                i2min,i2max=0,self.nx-dx-1
            else:
                i1min,i1max=0,self.nx+dx-1
                i2min,i2max=-dx,self.nx-1
            if dy>=0:
                j1min,j1max=dy,self.ny-1
                j2min,j2max=0,self.ny-dy-1
            else:
                j1min,j1max=0,self.ny+dy-1
                j2min,j2max=-dy,self.ny-1
            zsel=np.where(self.znum==i,self.zmax,0.0)
            ztrk[j2min:j2max,i2min:i2max]+=zsel[j1min:j1max,i1min:i1max]

        return ztrk
                
    def __repr__(self):
       return "%s %dx%dx%d %s %.3f %d %s"%(self.fname,self.nx,self.ny,self.nz,self.nfd,self.texp,self.site_id,self.observer)

   
