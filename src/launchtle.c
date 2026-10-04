#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>
#include "sgdp4h.h"
#include "elset.h"
#include "sattools_util.h"
#include "satutl.h"

#define LIM 128
#define XKMPER  6378.135            /* Km per earth radii */
#define XMNPDA  1440.0              /* Minutes per day */
#define AE      1.0                 /* Earth radius in "chosen units". */
#define XKE     0.743669161e-1
#define CK2     5.413080e-4   /* (0.5 * XJ2 * AE * AE) */
extern double SGDP4_jd0;

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

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// DOY to MJD
double doy2mjd(int year,double doy)
{
  int month,k=2;
  double day;

  if (year%4==0 && year%400!=0)
    k=1;

  month=floor(9.0*(k+doy)/275.0+0.98);
  
  if (doy<32.0)
    month=1;

  day=doy-floor(275.0*month/9.0)+k*floor((month+9.0)/12.0)+30.0;

  return date2mjd(year,month,day);
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}

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

void usage(void)
{
  printf("launch tle c:i:t:T:I:d:\n\n");
  printf("-c  reference tle\n-i  reference satno\n-t  reference launch time\n-T  launch time\n-I  output satno\n-d  output desig\n");

  return;
}

int main(int argc,char *argv[])
{
  int arg=0,satno=0,satno1=84001;
  char tlefile[LIM];
  char nfd0[32],nfd1[32];
  char line1[70],line2[70];
  FILE *file;
  orbit_t orb;
  double mjd0,mjd1,h0,h1,dmjd,dh;
  char *env;
  char desig[]="14900A";
  
  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/classfd.tle",env);

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:t:T:I:d:"))!=-1) {
    switch (arg) {
      
    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 't':
      strcpy(nfd0,optarg);
      mjd0=nfd2mjd(nfd0);
      break;

    case 'T':
      strcpy(nfd1,optarg);
      mjd1=nfd2mjd(nfd1);
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'I':
      satno1=atoi(optarg);
      break;

    case 'd':
      strcpy(desig,optarg);
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

  // Open file
  file=fopen(tlefile,"rb");
  if (file==NULL) 
    fatal_error("File open failed for reading \"%s\"",tlefile);
  
  // Find elements
  read_elset(file,satno,&orb,NULL);
  fclose(file);

  // Difference between epoch and launch time
  dmjd=doy2mjd(orb.ep_year,orb.ep_day)-mjd0;
  
  // Difference in RAAN
  orb.ascn=RAD(modulo(gmst(mjd1)-gmst(mjd0)+DEG(orb.ascn),360.0));

  // New epoch
  mjd0=mjd1+dmjd;
  
  // Format epoch
  orb.ep_day=mjd2doy(mjd0,&orb.ep_year);

  // Update desig
  strcpy(orb.desig,desig);
  orb.satno=satno1;

  // Print output
  if (format_tle_report(&orb,line1,line2)!=0) return 1;
  printf("%s\n%s\n",line1,line2);

  return 0;
}
