#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "sgdp4h.h"
#include "satutl.h"
#include <getopt.h>

#define R2D     180.0/M_PI
#define D2R     M_PI/180.0
#define LIM 128
#define XKMPER  6378.135 
#define XMNPDA  1440.0              /* Minutes per day */
#define XKE     0.743669161e-1

int fgetline(FILE *file,char *s,int lim);
void mjd2date(double mjd,int *year,int *month,double *day);
double date2mjd(int year,int month,double day);
double mjd2doy(double mjd,int *yr);
void format_tle(orbit_t orb,char *line1,char *line2);
double modulo(double x,double y);

void convert(char *line)
{
  int no,so;
  int year,month,day,hour,min,sec;
  int foy,fom,fod,gap;
  float sma,incl,raan,ecc,argl,argp,a2m,mag,dm;
  double mjd,doy;
  orbit_t orb;
  char line1[70],line2[70],line0[70];
  double e,v;
  float perigee,apogee;

  dm=5.0*log10(1000.0/40000.0);

  // Read line
  sscanf(line,"%4d,%6d,%2d%2d%4d,%2d%2d%4d %2d%2d%2d,%4d,%f,%f,%f,%f,%f,%f,%f,%f",&no,&so,&fod,&fom,&foy,&day,&month,&year,&hour,&min,&sec,&gap,&sma,&incl,&raan,&ecc,&argl,&argp,&a2m,&mag);

  // Use date of first observed as a designation
  mjd=date2mjd(foy,fom,fod);
  doy=mjd2doy(mjd,&foy);
  if (foy>=2000)
    sprintf(orb.desig,"%02d%3.0fA",foy-2000,doy+500);
  else
    sprintf(orb.desig,"%02d%3.0fA",foy-1900,doy+500);

  // Find epoch
  mjd=date2mjd(year,month,day+hour/24.0+min/1440.0+sec/86400.0);
  orb.ep_day=mjd2doy(mjd,&orb.ep_year);
  
  // Set values
  orb.satno=no;
  orb.eqinc=incl*D2R;
  orb.argp=argp*D2R;
  orb.ascn=raan*D2R;
  orb.ecc=ecc;
  orb.bstar=0.0;
  orb.rev=XKE*pow(sma/XKMPER,-1.5)*XMNPDA/(2.0*M_PI);

  // Mean anomaly
  v=argl-argp;
  e=2*atan(sqrt((1.0-ecc)/(1.0+ecc))*tan(v*D2R/2.0));
  orb.mnan=modulo((e-ecc*sin(e)),2.0*M_PI);
  
  // Orbit
  perigee=sma*(1.0-ecc)-XKMPER;
  apogee=sma*(1.0+ecc)-XKMPER;

  sprintf(line0,"SO %6d                      %4.1f             %7.0fkm  x%7.0fkm",so,mag+dm,perigee,apogee);

  // Format line
  format_tle(orb,line1,line2);
  printf("%s\n%s\n%s\n",line0,line1,line2);

  return;
}

int main(int argc,char *argv[])
{
  FILE *file;
  char line[LIM],*fname;
  int arg=0;

  // Decode options
  while ((arg=getopt(argc,argv,"f:"))!=-1) {
    switch(arg) {

    case 'f':
      fname=optarg;
      break;

    default:
      break;
    }
  }

  file=fopen(fname,"r");
  while (fgetline(file,line,LIM)>0) 
    convert(line);
  fclose(file);

  return 0;
}

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
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
  sprintf(line1,"1 %05dU %-8s %2d%012.8f  .00000000  00000-0 %8s 0    0",orb.satno,orb.desig,orb.ep_year-2000,orb.ep_day,sbstar);
  sprintf(line2,"2 %05d %8.4f %8.4f %07.0f %8.4f %8.4f %11.8f    0",orb.satno,DEG(orb.eqinc),DEG(orb.ascn),1E7*orb.ecc,DEG(orb.argp),DEG(orb.mnan),orb.rev);

  // Compute checksums
  for (i=0,csum=0;i<strlen(line1);i++) {
    if (isdigit(line1[i]))
      csum+=line1[i]-'0';
    else if (line1[i]=='-')
      csum++;
  }
  sprintf(line1,"%s%d",line1,csum%10);
  for (i=0,csum=0;i<strlen(line2);i++) {
    if (isdigit(line2[i]))
      csum+=line2[i]-'0';
    else if (line2[i]=='-')
      csum++;
  }
  sprintf(line2,"%s%d",line2,csum%10);

  return;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}
