#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "sgdp4h.h"

#define LIM 128
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)
#define STDMAG 6.0
#define PASSMAX 1000

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

int ipass=0,npass;
struct map {
  long satno;
  double lat,lng;
  double mjd;
  float alt,timezone;
  float saltmin,altmin;
  int length,all;
  char nfd[LIM],tlefile[LIM],observer[32];
  char datadir[LIM],tledir[LIM];
  int site_id,plot;
} m;
struct point {
  double mjd;
  char nfd[LIM];
  xyz_t obspos,satpos,sunpos;
  double sra,sde,salt,sazi,azi,alt;
} *pt;
struct pass {
  int satno;
  double mjdrise,mjdmax,mjdset;
  char line[80],skymap[LIM],radio[80];
  float length;
} p[PASSMAX];

double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void mjd2date(double mjd,char *date,int length);
void usage();
void initialize_setup(void);
void get_site(int site_id);
void nfd_now(char *s);
void print_header(void);
void obspos_xyz(double,xyz_t *,xyz_t *);
void sunpos_xyz(double,xyz_t *,double *,double *);
double gmst(double);
double dgmst(double);
double modulo(double,double);
void equatorial2horizontal(double,double,double,double *,double *);
void horizontal2equatorial(double,double,double,double *,double *);
void compute_observer_and_solar_positions(void);

void compute_track(orbit_t orb)
{
  int i,flag=0;
  double jd;
  xyz_t satpos,satvel;
  double dx,dy,dz,rsun,rearth;
  double h,psun,pearth,psat;
  char state[10]="",state1[10]="";
  double r,ra,de,phase,azi,alt,alt1;
  float mag;
  double mjdrise,mjdmax,mjdset,mjdentry,mjdexit;
  int irise=-1,imax=-1,iset=-1,ientry=-1,iexit=-1;
  int i1=-1,i2=-1,i3=-1;
  int sunlit,sundown;

  for (i=0;i<m.length;i++) {
    // Sun altitude
    if (pt[i].salt<m.saltmin)
      sundown=1;
    else
      sundown=0;

    // Compute satellite position
    jd=pt[i].mjd+2400000.5;
    satpos_xyz(jd,&satpos,&satvel);

    // Sun position from satellite
    dx=-satpos.x+pt[i].sunpos.x;  
    dy=-satpos.y+pt[i].sunpos.y;
    dz=-satpos.z+pt[i].sunpos.z;
    
    // Distances
    rsun=sqrt(dx*dx+dy*dy+dz*dz);
    rearth=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
    h=rearth-XKMPER;

    // Angles
    psun=asin(696.0e3/rsun)*R2D;
    pearth=asin(6378.135/rearth)*R2D;
    psat=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;
    
    // Visibility state
    if (psat-pearth<-psun)
      strcpy(state,"eclipsed");
    else if (psat-pearth>-psun && psat-pearth<psun)
      strcpy(state,"umbra");
    else if (psat-pearth>psun)
      strcpy(state,"sunlit");
    if (strcmp(state,"eclipsed")!=0)
      sunlit=1;
    else
      sunlit=0;

    // Position differences
    dx=satpos.x-pt[i].obspos.x;  
    dy=satpos.y-pt[i].obspos.y;
    dz=satpos.z-pt[i].obspos.z;
    
    // Celestial position
    r=sqrt(dx*dx+dy*dy+dz*dz);
    ra=modulo(atan2(dy,dx)*R2D,360.0);
    de=asin(dz/r)*R2D;
    
    // Phase
    phase=acos(((pt[i].obspos.x-satpos.x)*(pt[i].sunpos.x-satpos.x)+(pt[i].obspos.y-satpos.y)*(pt[i].sunpos.y-satpos.y)+(pt[i].obspos.z-satpos.z)*(pt[i].sunpos.z-satpos.z))/(rsun*r))*R2D;

    // Magnitude
    if (strcmp(state,"sunlit")==0) 
      mag=STDMAG-15.0+5*log10(r)-2.5*log10(sin(phase*D2R)+(M_PI-phase*D2R)*cos(phase*D2R));
    else
      mag=15;

    // Horizontal position
    equatorial2horizontal(pt[i].mjd,ra,de,&azi,&alt);
    if (alt>0.0)
      pt[i].alt=alt;
    else
      pt[i].alt=0.0;
    pt[i].azi=modulo(azi+180.0,360.0);

    // Find all passes
    if (m.all==1) {
      if (i==0) { 
	if (alt>=m.altmin && flag==0) {
	  irise=i;
	  flag=1;
	} 
	if (imax==-1 && flag==1) {
	  imax=i;
	  flag=2;
	}
      } else if (i==m.length-1) {
	if (alt>=m.altmin && flag==0) {
	  irise=i;
	  flag=1;
	} 
	if (imax==-1 && flag==1) {
	  imax=i;
	  flag=2;
	}
	if (alt>=m.altmin && flag==2) {
	  iset=i;
	  flag=3;
	}
      } else if (i>0) {
	if (alt1<m.altmin && alt>=m.altmin && flag==0) {
	  irise=i;
	  flag=1;
	} else if (alt1>m.altmin && alt<=m.altmin && flag==2) {
	  iset=i;
	  flag=3;
	} else if (flag==1 && alt<alt1) {
	  imax=i-1;
	  flag=2;
	} 
      }
    } else {
      if (flag==0 && alt>=m.altmin && sundown==1 && sunlit==1) {
	irise=i;
	flag=1;
      } else if (flag==1 && (alt<alt1 || !(sundown==1 && sunlit==1))) {
	imax=i-1;
	flag=2;
      } else if (flag==2 && !(alt>=m.altmin && sundown==1 && sunlit==1)) {
	iset=i;
	flag=3;
      }
      /*
      if (sundown==1 && sunlit==1) {
	if (alt>=m.altmin && flag==0) {
	  irise=i;
	  flag=1;
	} else if (alt<=m.altmin && flag==2) {
	  iset=i;
	  flag=3;
	} else if (flag==1 && alt<alt1) {
	  imax=i-1;
	  flag=2;
	} 
      }
      */
    }

    if (flag==3) {
      i1=irise;
      i2=imax;
      i3=iset;

      p[ipass].satno=orb.satno;
      p[ipass].mjdrise=pt[i1].mjd;
      p[ipass].mjdmax=pt[i2].mjd;
      p[ipass].mjdset=pt[i3].mjd;
      p[ipass].length=86400.0*(pt[i3].mjd-pt[i1].mjd);

      sprintf(p[ipass].line,"%05d | %s  %3.0f/%2.0f | %.8s  %3.0f/%2.0f | %.8s  %3.0f/%2.0f | \n",orb.satno,
	      pt[i1].nfd,pt[i1].azi,pt[i1].alt,
	      pt[i2].nfd+11,pt[i2].azi,pt[i2].alt,
	      pt[i3].nfd+11,pt[i3].azi,pt[i3].alt);
      sprintf(p[ipass].radio,"%05d %s %s %4.0f\n",orb.satno,pt[i1].nfd,pt[i3].nfd,p[ipass].length);
      sprintf(p[ipass].skymap,"skymap -c %s -i %d -s %d -t %s -l %.0f",m.tlefile,orb.satno,m.site_id,pt[i1].nfd,p[ipass].length);

      flag=0;
      irise=-1;
      imax=-1;
      iset=-1;
      ientry=-1;
      iexit=-1;
      ipass++;
      npass++;
    }

    alt1=alt;
    strcpy(state1,state);
  }

  return;
}

int qsort_compare_mjdrise(const void *a,const void *b) 
{ 
  struct pass *pa=(struct pass *) a;
  struct pass *pb=(struct pass *) b;

  if (pa->mjdrise<pb->mjdrise)
    return -1;
  else if (pa->mjdrise>pb->mjdrise)
    return 1;
  else
    return 0;
} 

int main(int argc,char *argv[])
{
  int arg=0,radio=0;
  FILE *file;
  orbit_t orb;
  int imode,quiet=0;

  // Initialize setup
  initialize_setup();

  // Decode options
  while ((arg=getopt(argc,argv,"t:c:i:s:l:hS:A:aPqm:R"))!=-1) {
    switch (arg) {

    case 'R':
      radio=1;
      break;
      
    case 't':
      strcpy(m.nfd,optarg);
      m.mjd=nfd2mjd(m.nfd);
      break;

    case 'm':
      m.mjd=atof(optarg);
      mjd2date(m.mjd,m.nfd,0);
      break;

    case 'c':
      strcpy(m.tlefile,optarg);
      break;

    case 's':
      get_site(atoi(optarg));
      break;

    case 'i':
      m.satno=atoi(optarg);
      break;

    case 'l':
      m.length=atoi(optarg);
      if (strchr(optarg,'h')!=NULL)
	m.length*=3600;
      else if (strchr(optarg,'m')!=NULL)
	m.length*=60;
      else if (strchr(optarg,'d')!=NULL)
	m.length*=86400;
      break;

    case 'S':
      m.saltmin=atof(optarg);
      break;

    case 'A':
      m.altmin=atof(optarg);
      break;

    case 'a':
      m.all=1;
      break;

    case 'h':
      usage();
      return 0;
      break;
      
    case 'P':
      m.plot=1;
      break;

    case 'q':
      quiet=1;
      break;

    default:
      usage();
      return 0;
    }
  }

  // Allocate
  pt=(struct point *) malloc(sizeof(struct point)*m.length);

  // Compute observer positions
  compute_observer_and_solar_positions();

  // Reloop stderr
  freopen("/tmp/stderr.txt","w",stderr);

  // Open TLE file
  file=fopen(m.tlefile,"r");
  if (file==NULL)
    fatal_error("File open failed for reading %s\n",m.tlefile);

  // Loop over objects
  while (read_twoline(file,m.satno,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);

    if (imode==SGDP4_ERROR)
      continue;

    // Skip non LEO objects
    if (orb.rev>=10.0 || m.satno!=0)
      compute_track(orb);
  }
  npass=ipass;
  
  // Close
  fclose(file);
  fclose(stderr);

  // Sort passes
  qsort(p,npass,sizeof(struct pass),qsort_compare_mjdrise);

  // Output header
  if (quiet==0)
    print_header();

  // Print passes
  for (ipass=0;ipass<npass;ipass++) {
    if (radio==0)
      printf("%s",p[ipass].line);
    else if (radio==1)
      printf("%s",p[ipass].radio);
    if (m.plot==1)
      system(p[ipass].skymap);
  }

  // Free
  free(pt);
  
  return 0;
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
  printf("pass t:c:i:s:l:hS:A:aPqm:R\n\n");
  printf("t    date/time (yyyy-mm-ddThh:mm:ss.sss) [default: now]\n");
  printf("c    TLE catalog file [default: classfd.tle]\n");
  printf("i    satellite ID (NORAD) [default: all]\n");
  printf("s    site (COSPAR\n");
  printf("l    length [default: %d s]\n",m.length);
  printf("A    minimum satellite altitude [default: %.1f deg]\n",m.altmin);
  printf("S    maximum solar altitude [default: %.1f deg]\n",m.saltmin);
  printf("a    compute all passes [toggle; default: off]\n");
  printf("P    plot passes [toggle; default: off]\n");
  printf("m    MJD date/time\n");
  printf("q    no header [toggle; default: off]\n");
  printf("R    format output for radio passes [toggle; default: off]\n");
  printf("h    this help\n");
  
  return;
}

// Compute Date from Julian Day
void mjd2date(double mjd,char *date,int length)
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

  if (length==3)
    sprintf(date,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);
  else if (length==0)
    sprintf(date,"%04d-%02d-%02dT%02d:%02d:%02.0f",year,month,day,hour,min,sec);
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

// Initialize setup
void initialize_setup(void)
{
  char *env;

  // Default parameters
  m.satno=0;
  m.timezone=0.0;
  m.length=3600;
  nfd_now(m.nfd);
  m.mjd=nfd2mjd(m.nfd);
  m.saltmin=-6.0;
  m.altmin=10.0;
  m.all=0;
  m.plot=0;

  // Default settings
  strcpy(m.observer,"Unknown");
  m.site_id=0;

  // Get environment variables
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
  env=getenv("ST_TLEDIR");
  if (env!=NULL) {
    strcpy(m.tledir,env);
  } else {
    printf("ST_TLEDIR environment variable not found.\n");
  }
  sprintf(m.tlefile,"%s/classfd.tle",m.tledir);

  return;
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

// Print header
void print_header(void)
{
  printf("Observer: %s (%04d) [%+.4f, %+.4f, %.0fm]\n",m.observer,m.site_id,m.lat,m.lng,m.alt*1000.0);
  printf("Elements: %s\n",m.tlefile);
  printf("UT Date/Time: %s for %g h \n",m.nfd,m.length/3600.0);
  return;
}

// Observer position
void obspos_xyz(double mjd,xyz_t *pos,xyz_t *vel)
{
  double ff,gc,gs,theta,s,dtheta;

  s=sin(m.lat*D2R);
  ff=sqrt(1.0-FLAT*(2.0-FLAT)*s*s);
  gc=1.0/ff+m.alt/XKMPER;
  gs=(1.0-FLAT)*(1.0-FLAT)/ff+m.alt/XKMPER;

  theta=gmst(mjd)+m.lng;
  dtheta=dgmst(mjd)*D2R/86400;

  pos->x=gc*cos(m.lat*D2R)*cos(theta*D2R)*XKMPER;
  pos->y=gc*cos(m.lat*D2R)*sin(theta*D2R)*XKMPER; 
  pos->z=gs*sin(m.lat*D2R)*XKMPER;
  vel->x=-gc*cos(m.lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
  vel->y=gc*cos(m.lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
  vel->z=0.0;

  return;
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

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
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

// Convert horizontal into equatorial coordinates
void horizontal2equatorial(double mjd,double azi,double alt,double *ra,double *de)
{
  double h;

  h=atan2(sin(azi*D2R),cos(azi*D2R)*sin(m.lat*D2R)+tan(alt*D2R)*cos(m.lat*D2R))*R2D;
  *ra=modulo(gmst(mjd)+m.lng-h,360.0);
  *de=asin(sin(m.lat*D2R)*sin(alt*D2R)-cos(m.lat*D2R)*cos(alt*D2R)*cos(azi*D2R))*R2D;
  if (*ra<0.0)
    *ra+=360.0;

  return;
}

void compute_observer_and_solar_positions(void)
{
  int i;
  xyz_t obsvel;

  for (i=0;i<m.length;i++) {
    // Compute MJDs
    pt[i].mjd=m.mjd+(double) i/86400.0;
    mjd2date(pt[i].mjd,pt[i].nfd,0);

    // Observer position
    obspos_xyz(pt[i].mjd,&pt[i].obspos,&obsvel);

    // Solar position
    sunpos_xyz(pt[i].mjd,&pt[i].sunpos,&pt[i].sra,&pt[i].sde);
    equatorial2horizontal(pt[i].mjd,pt[i].sra,pt[i].sde,&pt[i].sazi,&pt[i].salt);
  }

  return;
}
