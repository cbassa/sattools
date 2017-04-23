#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "sgdp4h.h"

#define LIM 384
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPAU 149597879.691 // AU in km

struct map {
  int site_id;
  double mjd;
  float saltmin,alt;
  char nfd[LIM],observer[32],datadir[LIM];
  double lat,lng;
  float length;
} m;
void allnight(void);
void sunpos_xyz(double,xyz_t *,double *,double *);
double modulo(double,double);
void equatorial2horizontal(double,double,double,double *,double *);
void mjd2date(double mjd,char *date);
double gmst(double mjd);
void nfd_now(char *s);
double nfd2mjd(char *date);
void get_site(int site_id);
double date2mjd(int year,int month,double day);
void usage(void);

int main(int argc,char *argv[])
{
  int arg=0;
  char *env;

  // Default settings
  strcpy(m.observer,"Unknown");
  m.site_id=0;
  m.mjd=-1.0;
  m.saltmin=-6.0;
  m.alt=0.0;
  m.lat=0.0;
  m.lng=0.0;
  
  // Get default site
  env=getenv("ST_DATADIR");
  if (env!=NULL) {
    strcpy(m.datadir,env);
  } else {
    printf("ST_DATADIR environment variable not found.\n");
  }
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    get_site(atoi(env));
  } else {
    printf("ST_COSPAR environment variable not found.\n");
  }

  // Get current time
  nfd_now(m.nfd);
  m.mjd=nfd2mjd(m.nfd);
  
  // Decode options
  while ((arg=getopt(argc,argv,"t:s:S:"))!=-1) {
    switch (arg) {

    case 't':
      strcpy(m.nfd,optarg);
      m.mjd=nfd2mjd(m.nfd);
      break;

    case 'S':
      m.saltmin=atof(optarg);
      break;

    case 's':
      get_site(atoi(optarg));
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
  
  // Compute set/rise times of sun
  allnight();

  return 0;
}

// Usage
void usage(void)
{
  printf("allnight t:s:S:\n\n");
  printf("t    date/time (yyyy-mm-ddThh:mm:ss.sss) [default: now]\n");
  printf("S    Minimum sun altitude\n");
  printf("s    site (COSPAR)\n");

  return;
}

// Convert equatorial into horizontal coordinates
void equatorial2horizontal(double mjd,double ra,double de,double *azi,double *alt)
{
  double h;

  h=gmst(mjd)+m.lng-ra;
  
  *azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(m.lat*D2R)-tan(de*D2R)*cos(m.lat*D2R))*R2D,360.0);
  *alt=asin(sin(m.lat*D2R)*sin(de*D2R)+cos(m.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

  return;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Solar position
void sunpos_xyz(double mjd,xyz_t *pos,double *ra,double *de)
{
  double jd,t,l0,m,e,c,r;
  double n,s,ecl;

  jd=mjd+2400000.5;
  t=(jd-2451545.0)/36525.0;
  l0=modulo(280.46646+t*(36000.76983+t*0.0003032),360.0)*D2R;
  m=modulo(357.52911+t*(35999.05029-t*0.0001537),360.0)*D2R;
  e=0.016708634+t*(-0.000042037-t*0.0000001267);
  c=(1.914602+t*(-0.004817-t*0.000014))*sin(m)*D2R;
  c+=(0.019993-0.000101*t)*sin(2.0*m)*D2R;
  c+=0.000289*sin(3.0*m)*D2R;

  r=1.000001018*(1.0-e*e)/(1.0+e*cos(m+c));
  n=modulo(125.04-1934.136*t,360.0)*D2R;
  s=l0+c+(-0.00569-0.00478*sin(n))*D2R;
  ecl=(23.43929111+(-46.8150*t-0.00059*t*t+0.001813*t*t*t)/3600.0+0.00256*cos(n))*D2R;

  *ra=atan2(cos(ecl)*sin(s),cos(s))*R2D;
  *de=asin(sin(ecl)*sin(s))*R2D;

  pos->x=r*cos(*de*D2R)*cos(*ra*D2R)*XKMPAU;
  pos->y=r*cos(*de*D2R)*sin(*ra*D2R)*XKMPAU;
  pos->z=r*sin(*de*D2R)*XKMPAU;

  return;
}

void allnight(void)
{
  int flag;
  xyz_t sunpos;
  double ra,de,azi,alt,alt0;
  double mjd,mjdrise=-1.0,mjdset=-1.0;
  char nfd[32];

  // Find solar altitude at reference time
  sunpos_xyz(m.mjd,&sunpos,&ra,&de);
  equatorial2horizontal(m.mjd,ra,de,&azi,&alt);

  // Sun below limit, find rise, then set
  if (alt<m.saltmin) {
    for (flag=0,mjd=m.mjd;mjd<m.mjd+0.5;mjd+=1.0/86400) {
      sunpos_xyz(mjd,&sunpos,&ra,&de);
      equatorial2horizontal(mjd,ra,de,&azi,&alt);
      
      if (flag!=0) {
	if (alt>m.saltmin && alt0<=m.saltmin)
	  mjdrise=mjd;
      }
    
      if (flag==0)
	flag=1;
      
      alt0=alt;
    }
    for (flag=0,mjd=m.mjd-0.5;mjd<m.mjd;mjd+=1.0/86400) {
      sunpos_xyz(mjd,&sunpos,&ra,&de);
      equatorial2horizontal(mjd,ra,de,&azi,&alt);
      
      if (flag!=0) {
	if (alt<m.saltmin && alt0>=m.saltmin)
	  mjdset=mjd;
      }
    
      if (flag==0)
	flag=1;
      
      alt0=alt;
    }
    // Sun above limit, find set, and rise
  } else {
    for (flag=0,mjd=m.mjd;mjd<m.mjd+1.0;mjd+=1.0/86400) {
      sunpos_xyz(mjd,&sunpos,&ra,&de);
      equatorial2horizontal(mjd,ra,de,&azi,&alt);
      
      if (flag!=0) {
	if (alt>m.saltmin && alt0<=m.saltmin)
	  mjdrise=mjd;
	if (alt<m.saltmin && alt0>=m.saltmin)
	  mjdset=mjd;
      }
    
      if (flag==0)
	flag=1;
      
      alt0=alt;
    }
  }

  m.mjd=mjdset;
  mjd2date(m.mjd,m.nfd);
  mjd2date(mjdrise,nfd);
  printf("%s %s\n",m.nfd,nfd);
  
  return;
}

// Compute Date from Julian Day
void mjd2date(double mjd,char *date)
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
  x=3600.*fabs(x)+0.0001;
  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  hour=x;
  sec=floor(1000.0*sec)/1000.0;

  sprintf(date,"%04d-%02d-%02dT%02d:%02d:%02.0f",year,month,day,hour,min,sec);

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

// Get observing site
void get_site(int site_id)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64],filename[LIM];

  sprintf(filename,"%s/data/sites.txt",m.datadir);
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with site information not found!\n");
    return;
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip
    if (strstr(line,"#")!=NULL)
      continue;

    // Strip newline
    line[strlen(line)-1]='\0';

    // Read data
    sscanf(line,"%4d %2s %lf %lf %f",
	   &id,abbrev,&lat,&lng,&alt);
    strcpy(observer,line+38);

    // Change to km
    alt/=1000.0;
    
    if (id==site_id) {
      m.lat=lat;
      m.lng=lng;
      m.alt=alt;
      m.site_id=id;
      strcpy(m.observer,observer);
    }

  }
  fclose(file);

  return;
}
