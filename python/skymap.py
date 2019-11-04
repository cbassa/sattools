#!/usr/bin/env python3
import numpy as np
import ppgplot
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, FK5, SkyCoord
from astropy.io import fits
from astropy import wcs

class Stars:
    """Tycho2 catalog"""

    def __init__(self):
#        hdu = fits.open("/home/bassa/code/c/satellite/sattools/data/tyc2.fits")
        hdu = fits.open("/home/bassa/hip.fits")

        self.ra = hdu[1].data.field('RA')*u.deg
        self.dec = hdu[1].data.field('DEC')*u.deg
        self.mag = hdu[1].data.field('MAG_VT')
        hdu.close()
        
        self.p = SkyCoord(ra=self.ra, dec=self.dec, frame="icrs")
        

class Skymap:

    def __init__(self):
        self.az = 00.0
        self.alt = 90.0
        self.w = 120.0
        self.level = 1
        self.minmag = -2.0
        self.maxmag = 5.0
        self.minrad = 0.02
        self.maxrad = 2.0
        self.orientation = "horizontal"
        self.sunalt = -20.0
        self.t = Time.now()
        self.site = 4171
        self.lat = 52.8344
        self.lon = 6.3785
        self.elev = 10.0
        self.camera = "A1600/55"
        self.length = 60.0
        self.fw = 18.3/2
        self.fh = 13.9/2
        self.observer = "Cees Bassa"
        self.q = 0.0
        self.loc = EarthLocation(lat=self.lat*u.deg, lon=self.lon*u.deg, height=self.elev*u.m)

    def update(self):
        self.t = Time.now()
        p = AltAz(az=self.az*u.deg, alt=self.alt*u.deg, obstime=self.t, location=self.loc).transform_to(FK5(equinox=self.t))
        self.ra = p.ra
        self.dec = p.dec

    def renew(self):
        if self.w>120.0:
            self.w = 120.0

        if self.level>9:
            self.level = 9
        if self.level<1:
            self.level = 1
            
        w = [120.0, 90.0, 60.0, 30.0, 20.0, 10.0, 5.0, 2.0, 1.0]
        minmag = [-2.0, -1.5, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0, 3.0]
        maxmag = [5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0]

        m.w = w[m.level-1]
        m.minmag = minmag[m.level-1]
        m.maxmag = maxmag[m.level-1]
        
    def plot_horizontal_grid(self):
        sch, sls, sci = ppgplot.pgqch(), ppgplot.pgqls(), ppgplot.pgqci()

        ppgplot.pgsci(15)
        ppgplot.pgsls(2)

        # Altitudes
        for alt in np.arange(0.0, 90.0, 20.0):
            az = np.arange(0.0, 360.0)
            rx, ry = forward(self.az, self.alt, az, alt*np.ones_like(az))
            for i in range(len(rx)):
                if i==0:
                    ppgplot.pgmove(rx[i], ry[i])
                if (np.abs(rx[i])<=1.5*self.w) & (np.abs(ry[i])<=self.w):
                    ppgplot.pgdraw(rx[i], ry[i])
                else:
                    ppgplot.pgmove(rx[i], ry[i])

        # Azimuths
        for az in np.arange(0.0, 360.0, 30.0):
            alt = np.arange(0.0, 80.0)
            rx, ry = forward(self.az, self.alt, az*np.ones_like(alt), alt)
            for i in range(len(rx)):
                if i==0:
                    ppgplot.pgmove(rx[i], ry[i])
                if (np.abs(rx[i])<=1.5*self.w) & (np.abs(ry[i])<=self.w):
                    ppgplot.pgdraw(rx[i], ry[i])
                else:
                    ppgplot.pgmove(rx[i], ry[i])
        ppgplot.pgsci(1)
        ppgplot.pgsls(1)
                    
        # Horizon
        az = np.arange(0.0, 360.0)
        rx, ry = forward(self.az, self.alt, az, np.zeros_like(az))
        for i in range(len(rx)):
            if i==0:
                ppgplot.pgmove(rx[i], ry[i])
            if (np.abs(rx[i])<=1.5*self.w) & (np.abs(ry[i])<=self.w):
                ppgplot.pgdraw(rx[i], ry[i])
            else:
                ppgplot.pgmove(rx[i], ry[i])    
            
def forward(l0, b0, l, b):
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---STG", "DEC--STG"]
    w.wcs.cd = [[1.0, 0.0], [0.0, 1.0]]
    w.wcs.crval = [l0, b0]
    w.wcs.crpix = [0.0, 0.0]
    
    return w.wcs_world2pix(np.stack((l, b), axis=-1), 1).T

def reverse(l0, b0, x, y):
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---STG", "DEC--STG"]
    w.wcs.cd = [[1.0, 0.0], [0.0, 1.0]]
    w.wcs.crval = [l0, b0]
    w.wcs.crpix = [0.0, 0.0]

    if type(x) is np.ndarray:
        l, b = w.wcs_pix2world(np.stack((x, y), axis=-1), 1).T
    else:
        l, b = w.wcs_pix2world([[x, y]], 1)[0]
    
    return l, b

if __name__ == "__main__":
    # Load stars
    s = Stars()

    # Load the master class
    m = Skymap()

    # Initialize plot
    ppgplot.pgopen("/xs")
    ppgplot.pgslw(2)
    ppgplot.pgpap(0.0, 0.75)

    # For ever loop
    redraw = True
    while True:
        # Redraw
        if redraw==True:
            # Update
            m.update()
            
            # Initialize window
            ppgplot.pgscr(0, 0., 0., 0.)
            ppgplot.pgeras()
            ppgplot.pgsvp(0.01, 0.99, 0.01, 0.99)
            ppgplot.pgwnad(-1.5*m.w, 1.5*m.w, -m.w, m.w)
        
            # Set background depending on solar altitude
            if m.sunalt>0.0:
                ppgplot.pgscr(0, 0.0, 0.0, 0.4)
            elif m.sunalt>-6.0:
                ppgplot.pgscr(0, 0.0, 0.0, 0.3)
            elif m.sunalt>-12.0:
                ppgplot.pgscr(0, 0.0, 0.0, 0.2)
            elif m.sunalt>-18.0:
                ppgplot.pgscr(0, 0.0, 0.0, 0.1)
            else:
                ppgplot.pgscr(0, 0.0, 0.0, 0.0)
            ppgplot.pgsci(0)
                
            # Plot box
            ppgplot.pgrect(-1.5*m.w, 1.5*m.w, -m.w, m.w)
            ppgplot.pgsci(1)
            ppgplot.pgsch(1.0)
            ppgplot.pgbox("BC", 0., 0, "BC", 0., 0)
            ppgplot.pgpt1(0., 0., 2)

            # Plot field-of-view
            ppgplot.pgsfs(2)
            ppgplot.pgrect(-m.fw, m.fw, -m.fh, m.fh)
            ppgplot.pgsfs(1)
        
            # Top left string
            ppgplot.pgsch(0.8)
            text = "%s UTC; %s (%04d) [%+.4f\\u\\(2218)\\d, %+.4f\\u\\(2218)\\d, %.0fm]" % (m.t.isot,
                                                                                            m.observer,
                                                                                            m.site,
                                                                                            m.lat,
                                                                                            m.lon,
                                                                                            m.elev)
            ppgplot.pgmtxt("T", 0.6, 0., 0., text)
        
            # Bottom string
            sra = m.ra.to_string(sep=":", unit="hourangle", pad=True, precision=2)
            sdec = m.dec.to_string(sep=":", unit="deg", pad=True, precision=1)
            text = "R: %s; D: %s; A: %.1f; E: %.1f; q: %.2f; S: %.1fx%.1f deg; " % (sra, sdec, m.az, m.alt, m.q, 3.0*m.w,2.0*m.w)
            text = text + "L: %d; O: %s; m < %.1f; C: %s; l: %.0f s" % (m.level, m.orientation, m.maxmag, m.camera, m.length)
            ppgplot.pgmtxt("B", 1.0, 0.0, 0.0, text)
            ppgplot.pgsch(1.0)

            # Plot stars
            c = s.mag<m.maxmag
            aa = s.p[c].transform_to(AltAz(location=m.loc, obstime=m.t))
            mag = s.mag[c]
            rad = (m.maxrad+(m.minrad-m.maxrad)*(mag-m.minmag)/(m.maxmag-m.minmag))*m.w/90.0
            rx, ry = forward(m.az, m.alt, aa.az.deg, aa.alt.deg)
            for i in range(len(rad)):
                ppgplot.pgsci(0)
                ppgplot.pgcirc(rx[i], ry[i], 1.3*rad[i])
                ppgplot.pgsci(1)
                ppgplot.pgcirc(rx[i], ry[i], rad[i])

            # Plot grid
            m.plot_horizontal_grid()
                
            # Reset redraw
            redraw = False
        
        # Get cursor
        x, y, char = ppgplot.pgcurs()

        print(x, y, char)
        
        # Quit
        if char==b"q":
            break

        # Reset
        if char==b"r":
            redraw = True

        # Recenter
        if char==b"c":
            m.az, m.alt = reverse(m.az, m.alt, x, y)
            redraw = True

        # Zoom in
        if ((char==b"+") | (char==b"=")) & (m.level<9):
            m.level+=1
            redraw = True

        # Zoom out
        if (char==b"-") & (m.level>1):
            m.level-=1
            redraw = True

        # Zenit
        if char==b"z":
            m.az, m.alt, m.level = 0.0, 90.0, 1
            redraw = True

        # South
        if char==b"s":
            m.az, m.alt, m.level = 180.0, 45.0, 3
            redraw = True

        # North
        if char==b"n":
            m.az, m.alt, m.level = 0.0, 45.0, 3
            redraw = True

        # East
        if char==b"e":
            m.az, m.alt, m.level = 90.0, 45.0, 3
            redraw = True

        # West
        if char==b"w":
            m.az, m.alt, m.level = 270.0, 45.0, 3
            redraw = True
            
            
        # Renew plot
        m.renew()

    # End
    ppgplot.pgend()
