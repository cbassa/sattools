#!/usr/bin/env python
import numpy as np
from astropy.time import Time
from astropy.io import fits
import time
import glob
import os

def read_pgm(fname):
    # Open file
    f=open(fname,"r")

    # Read lines
    line1=f.readline()
    line2=f.readline()
    line3=f.readline()
    line4=f.readline()

    # Read parameters
    nfd=line2.split(" ")[1]
    nx,ny=int(line3.split(" ")[0]),int(line3.split(" ")[1])

    # Read image
    z=np.fromfile(f,dtype='uint8',count=nx*ny).astype('float32')

    # Close file
    f.close()
    
    return z,nfd

# Create fits file
def create_fits_file(nx,ny,nz,path,iframe):
    # Allocate
    z=np.empty(nz*ny*nx,dtype='float32').reshape(nz,nx*ny)
    mjd=np.empty(nz,dtype='float64')
    dt=np.empty(nz,dtype='float32')

    # Loop over frames
    tstart=time.time()
    for i in xrange(nz):
        z[i],nfd=read_pgm(path+"/img%06d.pgm"%(iframe+i))

        # Format time
        t=Time(nfd,format='isot')
        if i==0:
            t0=t
        mjd[i]=t.mjd
        dt[i]=86400.0*(mjd[i]-mjd[0])
    print "Files read in %.3f s"%(time.time()-tstart)

    # Compute statistics
    tstart=time.time()
    zmax=np.max(z,axis=0)
    znum=np.argmax(z,axis=0)
    z1=np.sum(z,axis=0)-zmax
    z2=np.sum(z*z,axis=0)-zmax*zmax
    zavg=z1/float(nz-1)
    zstd=np.sqrt((z2-z1*zavg)/float(nz-2))
    print "Statistics in %.3f s"%(time.time()-tstart)

    # Reshape, reformat and flip
    tstart=time.time()
    zmax=np.flipud(zmax.astype('float32').reshape(ny,nx))
    znum=np.flipud(znum.astype('float32').reshape(ny,nx))
    zavg=np.flipud(zavg.astype('float32').reshape(ny,nx))
    zstd=np.flipud(zstd.astype('float32').reshape(ny,nx))
    z=np.array([zavg,zstd,zmax,znum])
    print "Reshape in %.3f s"%(time.time()-tstart)
    
    # Filename
    fname="%s.fits"%t0

    # Format header
    hdr=fits.Header()
    hdr['DATE-OBS']="%s"%t0
    hdr['MJD-OBS']=t0.mjd
    hdr['EXPTIME']=dt[-1]-dt[0]
    hdr['NFRAMES']=nz
    hdr['CRPIX1']=float(nx)/2.0
    hdr['CRPIX2']=float(ny)/2.0
    hdr['CRVAL1']=0.0
    hdr['CRVAL2']=0.0
    hdr['CD1_1']=1.0
    hdr['CD1_2']=0.0
    hdr['CD2_1']=0.0
    hdr['CD2_2']=1.0
    hdr['CTYPE1']="RA---TAN"
    hdr['CTYPE2']="DEC--TAN"
    hdr['CUNIT1']="deg"
    hdr['CUNIT2']="deg"
    hdr['CRRES1']=0.0
    hdr['CRRES2']=0.0
    hdr['EQUINIX']=2000.0
    hdr['RADECSYS']="ICRS"
    hdr['COSPAR']=4171
    hdr['OBSERVER']="Cees Bassa"
    for i in xrange(nz):
        hdr['DT%04d'%i]=dt[i]
    for i in xrange(10):
        hdr['DUMY%03d'%i]=0.0

    # Write fits file
    hdu=fits.PrimaryHDU(data=z,header=hdr)
    hdu.writeto(fname,clobber=True)

    return fname

# Main function
if __name__ == '__main__':
    # Settings
    nx=720
    ny=576
    nz=250
    path="/dev/shm/video0"

    # Start forever loop
    while True:
        # Find files
        files=sorted(glob.glob(path+"/img??????.pgm"))

        # Enough files
        if len(files)>nz:
            # Get first frame
            iframe=int(files[0].replace(path+"/img","").replace(".pgm",""))

            tstart=time.time()
            fname=create_fits_file(nx,ny,nz,path,iframe)
            tend=time.time()
            print "Created %s in %.2f"%(fname,tend-tstart)
            
            # Remove files
            for i in xrange(nz):
                os.remove(path+"/img%06d.pgm"%(iframe+i))
        else:
            time.sleep(1)
