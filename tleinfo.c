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
#define D2R M_PI/180.0
#define R2D 180.0/M_PI

extern double SGDP4_jd0;

void usage(void)
{
  printf("usage: tleinfo [-c TLEFILE] [-u] [|-1|-f] [ |-n|-d] [-i SATNO] [ |-a|-b] [-H] [-h]\n\n");

  printf("-c TLEFILE  The file containing orbital elements in the form of TLEs, 3 lines per object\n");
  printf("            default: ./bulk.tle\n");
  printf("-i SATNO    Filter only elements for objects with this satno\n");
  printf("-u          Show only one object (MODE0 only)\n");

  printf("\nSelect MODE:\n");
  printf("            MODE0: Show TLEs, object names or COSPAR designations\n");
  printf("-1          MODE1: Show list of elements (one line per object)\n");
  printf("-f          MODE2: Show human-readable parameters\n\n");

  printf("Select INFOTYPE\n");
  printf("MODE0:\n");
  printf("  default   Show the TLEs itself\n");
  printf("  -n        Show only the name of the objects\n");
  printf("  -d        Show only the COSPAR designation of the objects\n\n");
  printf("MODE1:\n");
  printf("  default   SATNO, YEAR, DOY, INCL, ASCN, ARGP, MA, ECC, MM, BSTAR\n");
  printf("  -a        SATNO, SEMI, PERIGEE, APOGEE, INCL, PERIOD, ECC\n");
  printf("  -b        SATNO, YEAR, DOY, INCL, ASCN, ARGP, MA, ECC, MM, floor(MJD), LNG_AT_MIDNIGHT\n\n");


  printf("-H          Show header (MODE1 only), default: disabled\n");
  printf("-h          Print usage\n");
}

double modulo(double x,double y);

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

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;
 
  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  //  if (c == '\n')
  //    s[i++] = c;
  s[i] = '\0';
  return i;
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

// Compute Date from Julian Day
void mjd2nfd(double mjd,char *nfd)
{
  double f,jd,dday;
  int z,alpha,a,b,c,d,e;
  int year,month,day,hour,min;
  float sec,x;

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

  dday=b-d-floor(30.6001*e)+f;
  if (e<14)
    month=e-1;
  else
    month=e-13;

  if (month>2)
    year=c-4716;
  else
    year=c-4715;

  day=(int) floor(dday);
  x=24.0*(dday-day);
  x=3600.*fabs(x);
  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  hour=x;
  sec=floor(1000.0*sec)/1000.0;

  sprintf(nfd,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);

  return;
}

float orbital_longitude_at_midnight(orbit_t orb,double mjd0)
{
  int rv,imode;
  double jd,tsince,mjd;
  kep_t K;

  imode=init_sgdp4(&orb);

  mjd=floor(mjd0);
  
  jd=mjd+2400000.5;
  tsince=1440.0*(jd-SGDP4_jd0);
  rv=sgdp4(tsince,1,&K);

  return modulo(K.theta*R2D,360.0);
}

int main(int argc,char *argv[])
{
  int arg=0,satno=0,header=0,oneline=0,no,name=0,desig=0;
  char tlefile[LIM];
  char line0[LIM],line1[LIM],line2[LIM],nfd[32];
  FILE *file;
  orbit_t orb;
  float aodp,perigee,apogee,period,lng;
  int info=0,unique=0;
  double mjd;
  char *env;

  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/bulk.tle",env);

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:aH1ftndbu"))!=-1) {
    switch (arg) {
      
    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 'u':
      unique=1;
      break;

    case '1':
      oneline=1;
      break;

    case 'f':
      oneline=2;
      break;

    case 'n':
      name=1;
      break;

    case 'd':
      desig=1;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'a':
      info=1;
      break;

    case 'b':
      info=2;
      break;

    case 'H':
      header=1;
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

  if (oneline==0) {
    // Open file
    file=fopen(tlefile,"rb");
    if (file==NULL) 
      fatal_error("File open failed for reading \"%s\"",tlefile);

    while (fgetline(file,line1,LIM)>0) {
      // Find TLE line
      if (line1[0]=='1') {
	    fgetline(file,line2,LIM);
	    sscanf(line1+2,"%d",&no);

	    if (satno==0 || satno==no) {
	      if (name==1 && desig==0)
	        printf("%s\n",line0);
	      else if (name==0 && desig==1)
	        printf("%.8s\n",line1+9);
	      else
	        printf("%s\n%s\n%s\n",line0,line1,line2);
	      if (unique==1)
	        break;
	    }
      }

      strcpy(line0,line1);
    }

    /*
    // Loop over file
    while (fgetline(file,line0,LIM)>0) {
      // Read data lines
      if (line0[0]!='1' || line0[0]!='2') {
	fgetline(file,line1,LIM);
	fgetline(file,line2,LIM);
	sscanf(line1+2,"%d",&no);

	if (satno==0 || satno==no) {
	  if (name==1 && desig==0) 
	    printf("%s\n",line0);
	  else if (name==0 && desig==1)
	    printf("%.8s\n",line1+9);
	  else
	    printf("%s\n%s\n%s\n",line0,line1,line2);
	}
      } else if (line0[0]=='1') {
	fgetline(file,line2,LIM);
	sscanf(line1+2,"%d",&no);

	if (satno==0 || satno==no)
	  printf("%s\n%s\n",line0,line2);
      }
    }
    */
    fclose(file);
  } else if (oneline==1) {
    // Open file
    file=fopen(tlefile,"rb");
    if (file==NULL) 
      fatal_error("File open failed for reading \"%s\"",tlefile);

    if (header==1) {
      if (info==0)
        printf("SATNO YEAR DOY     INCL    ASCN     ARGP     MA       ECC      MM   BSTAR\n");
      else if (info==1)
        printf("SATNO SEMI     PERIGEE  APOGEE    PERIOD  ECC\n");
      else if (info==2)
        printf("SATNO YEAR DOY     INCL    ASCN     ARGP     MA       ECC      MM   floor(MJD) LNG_AT_MIDNIGHT\n");
    }

    // Loop over file
    while (read_twoline(file,satno,&orb)==0) {
      orbit(orb,&aodp,&perigee,&apogee,&period);
      mjd=doy2mjd(orb.ep_year,orb.ep_day);
      mjd2nfd(mjd,nfd);
      if (info==0) printf("%05d %10.4lf %8.4f %8.4f %8.4f %8.4f %8.6f %8.5f %e\n",orb.satno,mjd,DEG(orb.eqinc),DEG(orb.ascn),DEG(orb.argp),DEG(orb.mnan),orb.ecc,orb.rev,orb.bstar);
      if (info==1) printf("%05d %6.0f x %6.0f x %6.2f %8.2f %8.6f %14.8lf\n",orb.satno,perigee,apogee,DEG(orb.eqinc),period,orb.ecc,mjd);
      if (info==2) {
	    lng=orbital_longitude_at_midnight(orb,mjd);
	    printf("%05d %10.4lf %8.4f %8.4f %8.4f %8.4f %8.6f %8.5f %10.4lf %8.4f\n",orb.satno,mjd,DEG(orb.eqinc),DEG(orb.ascn),DEG(orb.argp),DEG(orb.mnan),orb.ecc,orb.rev,floor(mjd),lng);
      }
    }
    fclose(file);
  } else if (oneline==2) {
    // Open file
    file=fopen(tlefile,"rb");
    if (file==NULL) 
      fatal_error("File open failed for reading \"%s\"",tlefile);

    // Loop over file
    while (read_twoline(file,satno,&orb)==0) 
      print_orb(&orb);

    fclose(file);
  }

  return 0;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}
