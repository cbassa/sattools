#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include "sgdp4h.h"
#include "satutl.h"
#include <getopt.h>

#define XKMPER  6378.135            /* Km per earth radii */
#define XMNPDA  1440.0              /* Minutes per day */
#define AE      1.0                 /* Earth radius in "chosen units". */
#define XKE     0.743669161e-1
#define CK2     5.413080e-4   /* (0.5 * XJ2 * AE * AE) */
#define R2D     180.0/M_PI
#define D2R     M_PI/180.0

extern double SGDP4_jd0;
void orbit(orbit_t orb,float *aodp,float *perigee,float *apogee,float *period);
void format_tle(orbit_t orb,char *line1,char *line2);
double mjd2doy(double mjd,int *yr);
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void mjd2date(double mjd,int *year,int *month,double *day);
double gmst(double mjd);
double modulo(double x,double y);

void usage(void)
{
  printf("faketle q:Q:i:I:w:t:m:n:d:\n\n");

  printf("-q   Perigee altitude (km)\n");
  printf("-Q   Apogee altitude (km)\n");
  printf("-I   Orbital inclination (deg)\n");
  printf("-n   RA of the ascending node (deg)\n");
  printf("-w   Argument of perigee (deg)\n");
  printf("-m   Mean anomaly (deg)\n");
  printf("-t   Epoch (YYYY-mm-ddThh:mm:ss)\n");
  printf("-i   Satellite number\n");
  printf("-d   Time offset from epoch (s)\n");

  return;
}

int main(int argc,char *argv[])
{
  orbit_t orb;
  float aodp,perigee,apogee,period,dt=0.0;
  char line1[70],line2[70];
  int arg=0;
  double mjd,dh;

  // Initialize
  orb.satno=99999;
  orb.eqinc=0.0;
  orb.ascn=0.0;
  orb.argp=0.0;
  orb.mnan=0.0;
  orb.bstar=0.5e-4;
  orb.ep_day=1.000;
  orb.ep_year=2013;

  // Decode options
  while ((arg=getopt(argc,argv,"q:Q:i:I:w:t:m:n:hd:"))!=-1) {
    switch(arg) {
    case 'q':
      perigee=atof(optarg);
      break;

    case 'Q':
      apogee=atof(optarg);
      break;

    case 'I':
      orb.eqinc=atof(optarg)*D2R;
      break;

    case 'n':
      orb.ascn=atof(optarg)*D2R;
      break;

    case 'i':
      orb.satno=atoi(optarg);
      break;

    case 'w':
      orb.argp=atof(optarg)*D2R;
      break;

    case 'm':
      orb.mnan=atof(optarg)*D2R;
      break;

    case 't':
      mjd=nfd2mjd(optarg);
      break;

    case 'd':
      dt=atof(optarg);
      break;

    case 'h':
      usage();
      return 0;

    default:
      usage();
      return 0;
    }
  }

  orb.ep_day=mjd2doy(mjd+dt/86400.0,&orb.ep_year);

  perigee+=XKMPER;
  apogee+=XKMPER;
  aodp=0.5*(perigee+apogee)/XKMPER;
  orb.ecc=0.5*(apogee-perigee)/(aodp*XKMPER);
  orb.rev=XKE*pow(aodp,-1.5)*XMNPDA/(2.0*M_PI);
  if (orb.rev<10)
    orb.bstar=0.0;
  orbit(orb,&aodp,&perigee,&apogee,&period);

  format_tle(orb,line1,line2);
  printf("%s\n%s\n",line1,line2);

  return 0;
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

// Format TLE
void format_tle(orbit_t orb,char *line1,char *line2)
{
  int i,csum;
  char sbstar[]=" 00000-0",bstar[13];

  // Format Bstar term
  if (fabs(orb.bstar)>1e-9) {
    sprintf(bstar,"%11.4e",10*orb.bstar);
    sbstar[0] = bstar[0];  sbstar[1] = bstar[1];  sbstar[2] = bstar[3];  sbstar[3] = bstar[4];
    sbstar[4] = bstar[5];  sbstar[5] = bstar[6];  sbstar[6] = bstar[8];  sbstar[7] = bstar[10];  sbstar[8] = '\0';
  }
  // Print lines
  sprintf(line1,"1 %05dU          %2d%012.8f  .00000000  00000-0 %8s 0    0",orb.satno,orb.ep_year-2000,orb.ep_day,sbstar);
  sprintf(line2,"2 %05d %8.4f %8.4f %07.0f %8.4f %8.4f %11.8f    0",orb.satno,DEG(orb.eqinc),DEG(orb.ascn),1E7*orb.ecc,DEG(orb.argp),DEG(orb.mnan),orb.rev);

  // Compute checksums
  for (i=0,csum=0;i<strlen(line1);i++) {
    if (isdigit(line1[i]))
      csum+=line1[i]-'0';
    else if (line1[i]=='-')
      csum++;
  }
  line1[strlen(line1) - 1] = '0' + (csum % 10);
  for (i=0,csum=0;i<strlen(line2);i++) {
    if (isdigit(line2[i]))
      csum+=line2[i]-'0';
    else if (line2[i]=='-')
      csum++;
  }
  line2[strlen(line2) - 1] = '0' + (csum % 10);

  return;
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

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}
