# Roadmap for `sattools` development
This document describes planned changes to the `sattools` repository.

## Background
The `sattools` began as a set of tools to help with tracking satellites using video and photographic cameras. The choice of `C` was one of familiarity with the code. Where possible, available libraries such as `wcslib`, `gsl`, `qfits` and `pgplot` were used. Some of these libraries are hard to build, raising the threshold for using `sattools`. 

With `python` becoming more and more popular, it makes sense to port some of the `sattools` functionality over to `python`, and use powerful libraries such as `numpy`, `scipy`, `astropy`, `matplotlib` and `opencv` for the core calculations of these tools. Porting the tools to another language has the added bonus that it is easier to solve mistakes made in designing the `C` code and improve functionality.

## Planned changes

The applications of the `sattools` repository can be divided into three categories; 
1. general purpose tools
1. tools for video observing
1. tools for photographic observing

Tools for radio observing have already been developed in their own repository `strf`: https://github.com/cbassa/strf. Porting tools for video observing is currently underway in the `stvid` repository: https://github.com/cbassa/stvid. The tools for photographic observing will live under a future repository, likely to be called `stphot`. An experimental project to use all-sky cameras for satellite observing is being developed in the `asm` repository: https://github.com/cbassa/asm (this could be renamed to `stasm` or `stas`; better suggestions are welcome).

Tools | `sattools` | `stvid` | `stphot` | Purpose | Future
------|------------|---------|----------|---------|--------
`addwcs` | | x | x | Add/fits a WCS to FITS file | Changed to a `python` function
`allnight` | x | | | Compute sunrise/sunset times | Stand alone tool
`angular` | x | | | Compute angular distances | Stand alone tool
`calibrate` | | x | x | Manually calibrate WCS | Stand alone tool
`confirm` | | | | Confirm satellite identifications | Deprecated, replaced with `stvid` functionality
`cvs2tle` | | | | | Deprecated
`dec2sex` | x | | | Convert decimal to sexagesimal | Stand alone tool
`deproject` | | | | | Deprecated
`detect` | | | | Detect satellite tracks | Deprecated, replaced with `stvid` functionality
`fakeiod` | x | | | Generate a fake IOD observation | Stand alone tool
`faketle` | x | | | Generate a fake TLE | Stand alone tool
`fitsheader` | x | | | Read a FITS header | Stand alone tool
`imgstat` | | | | | Deprecated, replaced with `stvid` functionality
`jpg2fits` | | | x | Convert a JPG image to FITS | Stand alone tool
`launchtle` | x | | | Adjust a TLE to a new launch time | Stand alone tool
`measure` | | | x | Measure satellite positions in photograhic observations | Stand alone tool
`mvtle` | | | | Move a TLE | Deprecated
`normal` | | | | Compute the normal of an orbit | Deprecated
`pass` | x | | | Compute satellite passes | Stand alone tool
`pgm2fits` | | | | Convert individual video frames to FITS | Deprecated, replaced with `stvid` functionality
`plotfits` | | | Manually calibrate WCS | Deprecated, replaced with `stvid` functionality
`posmatch` | | | | Match satellite positions | Deprecated
`posvel` | | | | Compute satellite position and velocity | Deprecated
`propagate` | x | | | Propagate a TLE to a new epoch | Stand alone tool

