#!/usr/bin/env python
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fname="2018-02-12T17:39:55.270.fits"

# Read fits
hdu=fits.open(fname)

# Read data
data=hdu[0].data
zavg,zstd,zmax,znum=data

# Read header
ny,nx=zavg.shape
nz=hdu[0].header['NFRAMES']
dt=[hdu[0].header['DT%04d'%i] for i in range(nz)]

# Sigma frame
zsig=(zmax-zavg)/(zstd+1e-9)

# Position
xm,ym=np.meshgrid(np.arange(nx),np.arange(ny))

# Selection width
w=10
sigma=5.0

# Selection
c=(zsig>sigma) & (xm>w) & (xm<nx-w) & (ym>w) & (ym<ny-w)

x=np.ravel(xm[c])
y=np.ravel(ym[c])
inum=np.ravel(znum[c]).astype('int')
sig=np.ravel(zsig[c])
t=np.array([dt[i] for i in inum])

# Load predictions
d=np.loadtxt(fname+".id",usecols=(1,2,3,4,5,6))

x0=d[:,0]
y0=d[:,1]
t0=np.zeros_like(x0)
x1=d[:,2]
y1=d[:,3]
t1=d[:,4]

for i in range(len(x0)):
    xr=x0[i]+(x1[i]-x0[i])*(t-t0[i])/(t1[i]-t0[i])
    yr=y0[i]+(y1[i]-y0[i])*(t-t0[i])/(t1[i]-t0[i])
    r=np.sqrt((x-xr)**2+(y-yr)**2)
    c=r<10.0
    print(i,np.sum(c))
    plt.figure()
    plt.plot([x0[i],x1[i]],[y0[i],y1[i]],'r')
    plt.scatter(x[c],y[c],c=t[c],s=10.0*(sig[c]-sigma))
    plt.scatter(x,y,c=t,s=1)
    plt.show()

#fig=plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#ax.scatter(x,y,t,c=t,s=10*(sig-5))

#for i in range(len(x0)):
#    ax.plot([x0[i],x1[i]],[y0[i],y1[i]],[t0[i],t1[i]],'k')
#plt.show()
