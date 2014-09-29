#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cel.h"
#include "sgdp4h.h"
#include <getopt.h>

#define LIM 80
#define NMAX 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

struct site {
  int id;
  double lng,lat;
  float alt;
  char observer[64];
};
struct site get_site(int site_id);
int fgetline(FILE *file,char *s,int lim);
double modulo(double x,double y);
double gmst(double mjd);
double dgmst(double mjd);
double date2mjd(int year,int month,double day);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel);
void mjd2date(double mjd,char *date);
double nfd2mjd(char *date);
void dec2sex(double x,char *s,int type);
double doy2mjd(int year,double doy);

void compute_position(double mjd,xyz_t satpos,struct site s,int satno,char *desig) 
{
  char sra[15],sde[15],nfd[32];
  xyz_t obspos,obsvel;
  double dx,dy,dz,mjd0=51544.5;
  double ra,de,ra0,de0,r;

  // Compute positions
  obspos_xyz(mjd,s.lng,s.lat,s.alt,&obspos,&obsvel);
  
  // compute difference vector
  dx=satpos.x-obspos.x;  
  dy=satpos.y-obspos.y;
  dz=satpos.z-obspos.z;
  
  // Celestial position
  r=sqrt(dx*dx+dy*dy+dz*dz);
  ra=modulo(atan2(dy,dx)*R2D,360.0);
  de=asin(dz/r)*R2D;
  
  // Precess position
  precess(mjd,ra,de,mjd0,&ra0,&de0);

  // Angle format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
  dec2sex(ra0/15.0,sra,0);
  dec2sex(de0,sde,1);

  // Get date
  mjd2date(mjd,nfd);

  // Format output
  printf("%05d %.2s %-6s %04d G %s 17 25 %s%s 37 S\n",satno,desig,desig+2,s.id,nfd,sra,sde);

  return;
}

int main(int argc,char *argv[])
{
  int arg=0,satno;
  struct site s;
  double mjd=0;
  char nfd[32],tlefile[LIM],*fname,line[LIM];
  int i,imode;
  FILE *file;
  orbit_t orb;
  xyz_t satpos,satvel;
  char *env;
  int usefile=0,usepos=0,status;

  // Get site
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    s=get_site(atoi(env));
  } else {
    printf("ST_COSPAR environment variable not found.\n");
  }

  // Decode options
  while ((arg=getopt(argc,argv,"t:c:i:s:f:p:"))!=-1) {
    switch (arg) {
      
    case 't':
      strcpy(nfd,optarg);
      mjd=nfd2mjd(nfd);
      break;

    case 'm':
      mjd=atof(optarg);
      break;

    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 's':
      s=get_site(atoi(optarg));
      break;

    case 'f':
      fname=optarg;
      usefile=1;
      break;

    case 'p':
      fname=optarg;
      usepos=1;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    default:
      return 0;
    }
  }

  // Get start mjd for finding elset
  if (usefile==1 && usepos==0) {
    file=fopen(fname,"r");
    fgetline(file,line,LIM);
    status=sscanf(line,"%lf",&mjd);
    fclose(file);
  }

  // Open catalog
  if (usepos==0) {
    file=fopen(tlefile,"r");
    if (file==NULL) 
      fatal_error("Failed to open %s\n",tlefile);
    
    // Read TLE
    do {
      status=read_twoline(file,satno,&orb);
    } while (doy2mjd(orb.ep_year,orb.ep_day)<mjd && status!=-1);
    fclose(file);

    // Check for match
    if (orb.satno!=satno) {
      //    fprintf(stderr,"object %d not found in %s\n",p.satno,filename);
      return;
    }
    
    // Initialize
    imode=init_sgdp4(&orb);
    if (imode==SGDP4_ERROR) {
      fprintf(stderr,"Error initializing SGDP4\n");
      exit(0);
    }
    
    // Compute
    if (usefile==0) {
      satpos_xyz(mjd+2400000.5,&satpos,&satvel);
      compute_position(mjd,satpos,s,orb.satno,orb.desig);
    } else {
      file=fopen(fname,"r");
      while (fgetline(file,line,LIM)>0) {
	status=sscanf(line,"%lf",&mjd);
	satpos_xyz(mjd+2400000.5,&satpos,&satvel);
	compute_position(mjd,satpos,s,orb.satno,orb.desig);
      }
      fclose(file);
    }
  } else {
    file=fopen(fname,"r");
    while (fgetline(file,line,LIM)>0) {
      status=sscanf(line,"%lf %lf %lf %lf",&mjd,&satpos.x,&satpos.y,&satpos.z);
      compute_position(mjd,satpos,s,99999,"14999A");
    }
    fclose(file);
  }
  return 0;
}

// Get observing site
struct site get_site(int site_id)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64];
  struct site s;
  char *env,filename[LIM];

  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/sites.txt",env);
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
    
    // Copy site
    if (id==site_id) {
      s.lat=lat;
      s.lng=lng;
      s.alt=alt;
      s.id=id;
      strcpy(s.observer,observer);
    }

  }
  fclose(file);

  return s;
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

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
}

// Observer position
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel)
{
  double ff,gc,gs,theta,s,dtheta;

  s=sin(lat*D2R);
  ff=sqrt(1.0-FLAT*(2.0-FLAT)*s*s);
  gc=1.0/ff+alt/XKMPER;
  gs=(1.0-FLAT)*(1.0-FLAT)/ff+alt/XKMPER;

  theta=gmst(mjd)+lng;
  dtheta=dgmst(mjd)*D2R/86400;

  pos->x=gc*cos(lat*D2R)*cos(theta*D2R)*XKMPER;
  pos->y=gc*cos(lat*D2R)*sin(theta*D2R)*XKMPER; 
  pos->z=gs*sin(lat*D2R)*XKMPER;
  vel->x=-gc*cos(lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
  vel->y=gc*cos(lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
  vel->z=0.0;

  return;
}

// Precess a celestial position
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de)
{
  double t0,t;
  double zeta,z,theta;
  double a,b,c;

  // Angles in radians
  ra0*=D2R;
  de0*=D2R;

  // Time in centuries
  t0=(mjd0-51544.5)/36525.0;
  t=(mjd-mjd0)/36525.0;

  // Precession angles
  zeta=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  zeta+=(0.30188-0.000344*t0)*t*t+0.017998*t*t*t;
  zeta*=D2R/3600.0;
  z=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  z+=(1.09468+0.000066*t0)*t*t+0.018203*t*t*t;
  z*=D2R/3600.0;
  theta=(2004.3109-0.85330*t0-0.000217*t0*t0)*t;
  theta+=-(0.42665+0.000217*t0)*t*t-0.041833*t*t*t;
  theta*=D2R/3600.0;
  
  a=cos(de0)*sin(ra0+zeta);
  b=cos(theta)*cos(de0)*cos(ra0+zeta)-sin(theta)*sin(de0);
  c=sin(theta)*cos(de0)*cos(ra0+zeta)+cos(theta)*sin(de0);

  *ra=(atan2(a,b)+z)*R2D;
  *de=asin(c)*R2D;

  if (*ra<360.0)
    *ra+=360.0;
  if (*ra>360.0)
    *ra-=360.0;

  return;
}

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\t')
    c=' ';
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
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

  sprintf(date,"%04d%02d%02d%02d%02d%05.0f",year,month,day,hour,min,sec*1000);

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

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int type)
{
  int i;
  double sec,deg,min,fmin;
  char sign;

  sign=(x<0 ? '-' : '+');
  x=60.*fabs(x);

  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;
  if (type==0)
    fmin=1000.0*(min-floor(min));
  else
    fmin=100.0*(min-floor(min));

  if (type==0)
    sprintf(s,"%02.0f%02.0f%03.0f",deg,floor(min),floor(fmin));
  else
    sprintf(s,"%c%02.0f%02.0f%02.0f",sign,deg,floor(min),floor(fmin));

  return;
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
