#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "sgdp4h.h"
#include "satutl.h"
#include <getopt.h>

#define LIM 128
#define XKE 0.07436680 // Guassian Gravitational Constant
#define XKMPER 6378.135
#define AE 1.0
#define XMNPDA 1440.0
#define R2D 180.0/M_PI
#define D2R M_PI/180.0

extern double SGDP4_jd0;


// Compute Date from Julian Day
void mjd2date(double mjd,int *year,int *month,double *day)
{
  double f,jd;
  int z,alpha,a,b,c,d,e;

  jd=mjd+2400000.5;
  jd+=0.5;

  z=floor(jd);
  f=fmod(jd,1.);

  if (z<2299161)
    a=z;
  else {
    alpha=floor((z-1867216.25)/36524.25);
    a=z+1+alpha-floor(alpha/4.);
  }
  b=a+1524;
  c=floor((b-122.1)/365.25);
  d=floor(365.25*c);
  e=floor((b-d)/30.6001);

  *day=b-d-floor(30.6001*e)+f;
  if (e<14)
    *month=e-1;
  else
    *month=e-13;

  if (*month>2)
    *year=c-4716;
  else
    *year=c-4715;

  return;
}

// MJD to DOY
double mjd2doy(double mjd,int *yr)
{
  int year,month,k=2;
  double day,doy;
  
  mjd2date(mjd,&year,&month,&day);

  if (year%4==0 && year%400!=0)
    k=1;

  doy=floor(275.0*month/9.0)-k*floor((month+9.0)/12.0)+day-30;

  *yr=year;

  return doy;
}


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
  if (year==1852 && month==10 && day<=4) b=0;

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

void usage(void)
{
  printf("propagate c:i:t:m:\n\nPropagates orbital elements to a new epoch using the SGP4/SDP4 model.\nDefault operation propagates classfd.tle to now,\n\n-c  input catalog\n-i  Satellite number\n-t  New epoch (YYYY-MM-DDTHH:MM:SS)\n-m  New epoch (MJD)\n");

  return;
}


int main(int argc,char *argv[])
{
  int imode,satno=0,arg,i;
  FILE *file;
  orbit_t orb;
  xyz_t r,v;
  char tlefile[LIM],nfd[32];
  char *env;
  double mjd;
  int length=3600,dl=60;

  // Get environment variable
  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/classfd.tle",env);

  // Set date
  nfd_now(nfd);
  mjd=nfd2mjd(nfd);

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:t:l:d:"))!=-1) {
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

    case 'l':
      length=atoi(optarg);
      break;

    case 'd':
      dl=atoi(optarg);
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

  // Reloop stderr
  freopen("/tmp/stderr.txt","w",stderr);

  // Open file
  file=fopen(tlefile,"r");
  while (read_twoline(file,satno,&orb)==0) {
    // Propagate
    imode=init_sgdp4(&orb);

    for (i=0;i<length;i+=dl) {
      satpos_xyz(mjd+2400000.5+(double) i/86400.0,&r,&v);
      printf("%f %f %f %f %f %f\n",r.x,r.y,r.z,v.x,v.y,v.z);
    }
  }
  fclose(file);

  return 0;
}
