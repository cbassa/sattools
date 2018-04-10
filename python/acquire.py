#!/usr/bin/env python
import sys
import numpy as np
import cv2
import time
import ctypes
import multiprocessing
from astropy.time import Time
from astropy.io import fits
from sun import get_sunriseset

# Capture images
def capture(buf,z1,t1,z2,t2,device,nx,ny,nz,tend):
    # Array flag
    first=True

    # Loop until reaching end time
    while float(time.time())<tend:
        # Get frames
        for i in range(nz):
            # Get frame
            ret,frame=device.read()
            
            # Skip lost frames
            if ret==True:
                # Store time
                t=float(time.time())
                
                # Convert image to grayscale
                gray=np.asarray(cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)).astype(np.uint8)
                
                # Store buffer        
                z=gray

                # Store results
                if first:
                    z1[i]=z
                    t1[i]=t
                else:
                    z2[i]=z
                    t2[i]=t

        # Assign buffer ready
        if first:
            buf.value=1
        else:
            buf.value=2
        
        # Swap flag
        first=not first

    print("Exiting capture")    
        
def compress(buf,z1,t1,z2,t2,nx,ny,nz,tend):
    # Flag to keep track of processed buffer
    process_buf=1

    # Start processing
    while True:
        # Wait for buffer to become available
        while buf.value!=process_buf:
            time.sleep(1.0)

        # Process first buffer
        if buf.value==1:
            t=t1
            z=z1.astype('float32')
            process_buf=2
        elif buf.value==2:
            t=t2
            z=z2.astype('float32')
            process_buf=1

        # Compute statistics
        zmax=np.max(z,axis=0)
        znum=np.argmax(z,axis=0)
        zs1=np.sum(z,axis=0)-zmax
        zs2=np.sum(z*z,axis=0)-zmax*zmax
        zavg=zs1/float(nz-1)
        zstd=np.sqrt((zs2-zs1*zavg)/float(nz-2))

        # Convert to float and flip
        zmax=np.flipud(zmax.astype('float32'))
        znum=np.flipud(znum.astype('float32'))
        zavg=np.flipud(zavg.astype('float32'))
        zstd=np.flipud(zstd.astype('float32'))

        # Format time
        nfd="%s.%03d"%(time.strftime("%Y-%m-%dT%T", time.gmtime(t[0])),int((t[0]-np.floor(t[0]))*1000))
        t0=Time(nfd,format='isot')
        dt=t-t[0]

        # Generate fits
        fname="%s.fits"%nfd

        # Format header
        hdr=fits.Header()
        hdr['DATE-OBS']="%s"%nfd
        hdr['MJD-OBS']=t0.mjd
        hdr['EXPTIME']=dt[-1]-dt[0]
        hdr['NFRAMES']=nz
        hdr['CRPIX1']=float(nx)/2.0
        hdr['CRPIX2']=float(ny)/2.0
        hdr['CRVAL1']=0.0
        hdr['CRVAL2']=0.0
        hdr['CD1_1']=1.0/3600.0
        hdr['CD1_2']=0.0
        hdr['CD2_1']=0.0
        hdr['CD2_2']=1.0/3600.0
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
        for i in range(nz):
            hdr['DT%04d'%i]=dt[i]
        for i in range(10):
            hdr['DUMY%03d'%i]=0.0
        
        # Write fits file
        hdu=fits.PrimaryHDU(data=np.array([zavg,zstd,zmax,znum]),header=hdr)
        hdu.writeto(fname)
        print("Compressed",fname)

        # Exit on end of capture
        if t[-1]>tend:
            break

    # Exiting
    print("Exiting compress")
        
 # Main function
if __name__ == '__main__':
    # Current time
    tnow=float(time.time())
    tnow=1520573400.0
    tnow=1520510226.281413,
    
    # Get sunrise and sunset times
    trise,tset=get_sunriseset(tnow,52.8344,6.3785,10,-6.0)

    print(tnow,trise,tset)
    # Wait until sunset
    if (tset<trise) & (tnow<tset):
        print("Waiting for sunset at %s"%time.strftime("%FT%T",time.gmtime(tset)))

    sys.exit()


    # Settings
    devid=int(sys.argv[1])
    nx=720
    ny=576
    nz=250

    # Initialize device
    device=cv2.VideoCapture(devid)

    # Set properties
    device.set(3,nx)
    device.set(4,ny)

    # Initialize arrays
    z1base=multiprocessing.Array(ctypes.c_uint8,nx*ny*nz)
    z1=np.ctypeslib.as_array(z1base.get_obj()).reshape(nz,ny,nx)
    t1base=multiprocessing.Array(ctypes.c_double,nz)
    t1=np.ctypeslib.as_array(t1base.get_obj())
    z2base=multiprocessing.Array(ctypes.c_uint8,nx*ny*nz)
    z2=np.ctypeslib.as_array(z2base.get_obj()).reshape(nz,ny,nx)
    t2base=multiprocessing.Array(ctypes.c_double,nz)
    t2=np.ctypeslib.as_array(t2base.get_obj())
    buf=multiprocessing.Value('i',0)

    
    tend=float(time.time())+31.0
    
    # Set processes
    pcapture=multiprocessing.Process(target=capture,args=(buf,z1,t1,z2,t2,device,nx,ny,nz,tend))
    pcompress=multiprocessing.Process(target=compress,args=(buf,z1,t1,z2,t2,nx,ny,nz,tend))

    # Start
    pcapture.start()
    pcompress.start()

    # End
    try:
        pcapture.join()
        pcompress.join()
    except KeyboardInterrupt:
        pcapture.terminate()
        pcompress.terminate()
                        
    # Release device
    device.release()


