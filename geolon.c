#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "sgdp4h.h"

#define LIM 128
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;
struct map {
  long satno;
  double mjd;
  char nfd[LIM],tlefile[LIM],observer[32];
  char datadir[LIM],tledir[LIM];
} m;


// Compute Julian Day from Date
double date2mjd(int year,int month,double day)
{
  int a,b;
  double jd;

  if (month<3) {
    year--;
    month+=12;
  }

  a=floor(year/100.);
  b=2.-a+floor(a/4.);

  if (year<1582) b=0;
  if (year==1582 && month<10) b=0;
  if (year==1582 && month==10 && day<=4) b=0;

  jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

  return jd-2400000.5;
}

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min,sec;
  double mjd,dday;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%02d",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;

  mjd=date2mjd(year,month,dday);

  return mjd;
}

void usage()
{
  printf("Usage goes here.\n");
  
  return;
}


// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}

// Present nfd
void nfd_now(char *s)
{
  time_t rawtime;
  struct tm *ptm;

  // Get UTC time
  time(&rawtime);
  ptm=gmtime(&rawtime);
    
  sprintf(s,"%04d-%02d-%02dT%02d:%02d:%02d",ptm->tm_year+1900,ptm->tm_mon+1,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
  
  return;
}

// Compute longitude
void compute_longitude(char *tlefile,long satno,double mjd)
{
  FILE *fp=NULL;
  orbit_t orb;
  xyz_t satpos,satvel;
  double jd,h,l,b,r;
  long imode;
  
  // Open TLE file
  fp = fopen(tlefile, "rb");
  if(fp == NULL) 
    fatal_error("File open failed for reading \"%s\"", tlefile);
  
  // Loop over elements
  while(read_twoline(fp, satno, &orb) == 0) {
    Isat = orb.satno;
    imode = init_sgdp4(&orb);
    if(imode == SGDP4_ERROR) continue;

    // Skip objects with mean motions outside of 0.8 to 1.2 revs/day
    if (orb.rev<0.8 || orb.rev>1.2)
      continue;
    
    // Get Julian Date
    jd=mjd+2400000.5;
    
    // Get positions
    satpos_xyz(jd,&satpos,&satvel);

    // Greenwich Mean Sidereal time
    h=gmst(mjd);

    // Celestial position
    r=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
    l=atan2(satpos.y,satpos.x)*R2D;
    l=modulo(l-h,360.0);
    b=asin(satpos.z/r)*R2D;
    if (l>180.0) 
      l-=360.0;
    if (l<-180.0) 
      l+=360.0;

    printf("%05d %10s %8.3lf %8.3lf %6.0lf\n",orb.satno,orb.desig,l,b,r-XKMPER);
  }
  fclose(fp);

  return;
}

int main(int argc,char *argv[])
{
  int arg=0;
  long satno=0;
  double mjd;
  char nfd[LIM],tlefile[LIM];

  nfd_now(nfd);
  mjd=nfd2mjd(nfd);

  
  // Decode options
  if (argc>1) {
    while ((arg=getopt(argc,argv,"t:c:i:h"))!=-1) {
      switch (arg) {
	
      case 't':
	strcpy(nfd,optarg);
	mjd=nfd2mjd(nfd);
	break;
	
      case 'c':
	strcpy(tlefile,optarg);
	break;
	
      case 'i':
	satno=atoi(optarg);
	break;
	
      case 'h':
	usage();
	return 0;
	break;
	
      default:
	usage();
	return 0;
      }
    }
  } else {
    usage();
    return 0;
  }    

  // Compute longitudes of satellites
  compute_longitude(tlefile,satno,mjd);

  return 0;
}
