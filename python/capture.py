#!/usr/bin/env python

import numpy as np
import cv2
import time
import subprocess
import os

# Make directory
if not os.path.exists("/dev/shm/video0"):
    os.mkdir("/dev/shm/video0")

# Get sunset and sunrise times (sattools allnight program)
tsolar=subprocess.check_output(['allnight']).replace("\n","").split(" ")

# Convert to time structs
tset=time.strptime(tsolar[0]+" UTC", "%Y-%m-%dT%H:%M:%S %Z")
trise=time.strptime(tsolar[1]+" UTC", "%Y-%m-%dT%H:%M:%S %Z")
tnow=time.gmtime()
dtset=time.mktime(tset)-time.mktime(tnow)
dtrise=time.mktime(trise)-time.mktime(tnow)

# Wait for sunset
#if dtset>0:
#    print time.strftime("%FT%T",time.gmtime())+" waiting for sunset (%ds)"%dtset
#    time.sleep(dtset)

# Settings
device=cv2.VideoCapture(0)
device.set(3,720)
device.set(4,576)

# Set counter
iframe=1

# Start capture
print time.strftime("%FT%T",time.gmtime())+" start capture"
while dtrise>0:
    # Get frame
    ret,frame=device.read()

    # Skip lost frames
    if ret==True:
        # Get time
        t=float(time.time())

        # Format time
        nfd="%s.%03d"%(time.strftime("%Y-%m-%dT%T", time.gmtime(t)),int((t-np.floor(t))*1000))
                
        # Get Size
        ny,nx=frame.shape[0],frame.shape[1]

        # Convert image to grayscale
        gray=np.asarray(cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)).astype(np.int8)

        # Open output file
        f=open("/dev/shm/video0/img%06d.pgm"%iframe,"w")
        f.write("P5\n# %s\n%d %d\n255\n"%(nfd,nx,ny))
        f.write(gray)
        f.close()

        # Log
        if iframe%250==1:
            print "%06d: %s %dx%d"%(iframe,nfd,nx,ny)
        
        # Increment
        iframe+=1
    
    # Refresh time
    tnow=time.gmtime()
    dtrise=time.mktime(trise)-time.mktime(tnow)

    
# End capture
print time.strftime("%FT%T",time.gmtime())+" end capture"
