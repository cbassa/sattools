#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "sgdp4h.h"
#include "satutl.h"

#define LIM 128
#define XKMPER  6378.135            /* Km per earth radii */
#define XMNPDA  1440.0              /* Minutes per day */
#define AE      1.0                 /* Earth radius in "chosen units". */
#define XKE     0.743669161e-1
#define CK2     5.413080e-4   /* (0.5 * XJ2 * AE * AE) */
extern double SGDP4_jd0;

void usage(void)
{
  return;
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

// DOY to MJD
double doy2mjd(int year,double doy)
{
  int month,k=2;
  double day;

  if (year%4==0 && year%400!=0)
    k=1;

  month=floor(9.0*(k+doy)/275.0+0.98);
  
  if (doy<32)
    month=1;

  day=doy-floor(275.0*month/9.0)+k*floor((month+9.0)/12.0)+30.0;

  return date2mjd(year,month,day);
}

void orbit(orbit_t orb,float *aodp,float *perigee,float *apogee,float *period)
{
  float xno,eo,xincl;
  float a1,betao2,betao,temp0,del1,a0,del0,xnodp;

  xno=orb.rev*2.0*M_PI/XMNPDA;
  eo=orb.ecc;
  xincl=orb.eqinc;

  a1 = pow(XKE / xno, 2.0/3.0);
  betao2 = 1.0 - eo * eo;
  betao = sqrt(betao2);
  temp0 = (1.5 * CK2) * cos(xincl)*cos(xincl) / (betao * betao2);
  del1 = temp0 / (a1 * a1);
  a0 = a1 * (1.0 - del1 * (1.0/3.0 + del1 * (1.0 + del1 * 134.0/81.0)));
  del0 = temp0 / (a0 * a0);
  xnodp = xno / (1.0 + del0);
  *aodp = (a0 / (1.0 - del0));
  *perigee = (*aodp * (1.0 - eo) - 1) * XKMPER;
  *apogee = (*aodp * (1.0 + eo) - 1) * XKMPER;
  *period = (TWOPI * 1440.0 / XMNPDA) / xnodp;
  *aodp=(*aodp-1)*XKMPER;

  return;
}

int main(int argc,char *argv[])
{
  int arg=0,satno=0,quiet=0;
  char tlefile[LIM];
  char line0[70],line1[70],line2[70];
  FILE *file;
  orbit_t orb;
  float aodp,perigee,apogee,period;
  int info=0;
  double mjd;
  char *env;

  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/classfd.tle",env);

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:aq"))!=-1) {
    switch (arg) {
      
    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'a':
      info=1;
      break;

    case 'q':
      quiet=1;
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

  if (info==0 && quiet==0)
    printf("SATNO YEAR DOY     INCL    ASCN     ARGP     MA       ECC      MM\n");
  if (info==1 && quiet==0)
    printf("SATNO SEMI     PERIGEE  APOGEE    PERIOD  ECC\n");

  // Loop over file
  while (read_twoline(file,satno,&orb)==0) {
    orbit(orb,&aodp,&perigee,&apogee,&period);
    mjd=doy2mjd(orb.ep_year,orb.ep_day);
    if (info==0) printf("%05d %10.4lf %8.4f %8.4f %8.4f %8.4f %8.6f %8.5f\n",orb.satno,mjd,DEG(orb.eqinc),DEG(orb.ascn),DEG(orb.argp),DEG(orb.mnan),orb.ecc,orb.rev);
    if (info==1) printf("%05d %9.2f %9.2f %9.2f %8.2f %8.6f %14.8lf\n",orb.satno,aodp,perigee,apogee,period,orb.ecc,mjd);
  }
  fclose(file);
 
  return 0;
}
