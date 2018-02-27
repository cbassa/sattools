import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import cv2
import sys

hdu=fits.open(sys.argv[1])

data=hdu[0].data
zavg,zstd,zmax,znum=data

zsig=(zmax-zavg)/(zstd+1e-9)

ny,nx=zsig.shape
tmp=np.zeros(nx*ny*3,dtype='uint8').reshape(ny,nx,3)
tmp[:,:,0]=zmax.astype('uint8')
tmp[:,:,1]=0.1*zmax.astype('uint8')
tmp[:,:,2]=0.1*zmax.astype('uint8')
img=cv2.cvtColor(tmp,cv2.COLOR_BGR2RGB)

mask=np.array(zsig>5.0,dtype='uint8')
tmp=np.zeros(nx*ny*3,dtype='uint8').reshape(ny,nx,3)
tmp[:,:,0]=255*mask
tmp[:,:,1]=255*mask
tmp[:,:,2]=255*mask

mat=cv2.cvtColor(tmp,cv2.COLOR_BGR2GRAY)
lines=cv2.HoughLines(mat,1,np.pi/180.0,200)
if lines is not None:
    for rho,theta in lines[0]:
        a = np.cos(theta)
        b = np.sin(theta)
        x0 = a*rho
        y0 = b*rho
        x1 = int(x0 + 1000*(-b))
        y1 = int(y0 + 1000*(a))
        x2 = int(x0 - 1000*(-b))
        y2 = int(y0 - 1000*(a))

        cv2.line(img,(x1,y1+10),(x2,y2+10),(0,0,255),2)

minLineLength = 100
maxLineGap = 10
lines = cv2.HoughLinesP(mat,1,np.pi/180,100,minLineLength,maxLineGap)
if lines is not None:
    for x1,y1,x2,y2 in lines[0]:
        cv2.line(img,(x1,y1),(x2,y2),(0,255,0),2)

    cv2.imwrite(sys.argv[1]+'.jpg',img)
