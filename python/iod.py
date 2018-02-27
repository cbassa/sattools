#!/usr/bin/env python

class iod:
    """IOD Observation"""

    def __init__(self,line):
        s=line.split()
        self.iodline=line
        self.norad=int(s[0])
        self.cospar="%s %s"%(s[1],s[2])
        self.site=int(s[3])
        self.conditions=s[4]
        self.datestring=s[5]


        
    def __repr__(self):
        return "%05d %s %04d %s %s"%(self.norad,self.cospar,self.site,self.conditions,self.datestring)
        

iodlines=["10529 77 112D   4553 G 20120511214446298 17 25 2034286+472232 37 R",
          "10529 77 112D   4553 G 20120511214452938 17 25 2038653+454274 37 S",
          "23862 96 029D   4553 G 20120511221706217 17 25 2121225+480700 37 S",
          "23862 96 029D   4553 G 20120511221726241 17 25 2122590+453398 37 R"]

obslist=[iod(line) for line in iodlines]
    
print(type(obslist),type(obslist[0]))
for obs in obslist:
    print(obs)
