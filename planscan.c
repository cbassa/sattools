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
#define XKE     0.743669161e-1
#define CK2     5.413080e-4   /* (0.5 * XJ2 * AE * AE) */
#define XMNPDA  1440.0              /* Minutes per day */
#define FLAT (1.0/298.257)
#define STDMAG 6.0
#define MMAX 1024


long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;
struct map {
  double lat,lng;
  float alt;
  char observer[32];
  int site_id;
} m;
void get_site(int site_id);
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void nfd_now(char *s);
void obspos_xyz(double mjd,xyz_t *pos,xyz_t *vel);
void sunpos_xyz(double mjd,xyz_t *pos,double *ra,double *de);
double gmst(double mjd);
double dgmst(double mjd);
double modulo(double x,double y);
void equatorial2horizontal(double mjd,double ra,double de,double *azi,double *alt);
void mjd2date(double mjd,char *date);
int properties(kep_t K,xyz_t obspos,xyz_t sunpos,float radius,float t,double *ra,double *de,double *r,float *mag);
void dec2sex(double x,char *s,int f,int len);

float semimajoraxis(orbit_t orb)
{
  float xno,eo,xincl;
  float a1,betao2,betao,temp0,del1,a0,del0,xnodp,aodp;

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
  aodp = (a0 / (1.0 - del0));
  aodp=(aodp-1)*XKMPER;

  return aodp;
}

void usage(void)
{
  printf("planscan -t <UT Date/time> -c <catalog> -s <site> -l <length> -i <NORAD>\n");
  printf("         -r <altitude> -A <elevation> -S <elevation>\n\n");
  printf("-t <UT Date/time>  UT Start date/time in yyyy-mm-ddThh:mm:ss [default: now]\n");
  printf("-c <catalog>       Input TLE catalog to use [default: classfd.tle]\n");
   printf("-s <site>          Site number from sites.txt [default: 4171]\n");
   printf("-l <length>        Search length from UT start in seconds [default: 86400 s]\n");
   printf("-i <NORAD>         NORAD number of satellite to select [default: 41334]\n");
   printf("-r <altitude>      Satellite altitude above surface in km [default: mean orbital altitude]\n");
   printf("-A <elevation>     Minimum satellite elevation in degrees [default: 10 degrees]\n");
   printf("-S <elevation>     Maximum solar elevation in degrees [default: -6 degrees\n");
   printf("-d <timestep>      Time step in seconds [default: 60s]\n");
   printf("-C                 Select on culmination instead of maximum brightness\n");
   printf("-h                 Shows this help\n");

  return;
}

int main(int argc,char *argv[])
{
  int arg=0,imode,i;
  int satno=-1;
  orbit_t orb;
  FILE *fp;
  char tlefile[256]="classfd.tle";
  char nfd[32];
  double mjd0,mjd,jd,tsince;
  int rv,withvel;
  kep_t K;
  double radius=-1;
  xyz_t obspos,obsvel,sunpos;
  double p,pmin,p0,p1,r,ra,de,azi,alt,sazi,salt,altmin=10.0,saltmin=-6.0,altmax,dp;
  float mag,mmin;
  int state,pstate,nstate;
  float t,length=86400.0,dt=60.0;
  char sra[16],sde[16],type[32];
  int opttype=0; // 0 magnitude; 1 elevation

  // Initialize
  nfd_now(nfd);
  mjd0=nfd2mjd(nfd);
  get_site(4171);
  
  // Decode options
  while ((arg=getopt(argc,argv,"t:c:i:s:l:hS:A:r:d:C"))!=-1) {
    switch (arg) {
      
    case 't':
      strcpy(nfd,optarg);
      mjd0=nfd2mjd(nfd);
      break;

    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 's':
      get_site(atoi(optarg));
      break;

    case 'C':
      opttype=1;
      break;
      
    case 'i':
      satno=atoi(optarg);
      break;

    case 'r':
      radius=atof(optarg);
      break;
      
    case 'l':
      length=atoi(optarg);
      if (strchr(optarg,'h')!=NULL)
	length*=3600;
      else if (strchr(optarg,'m')!=NULL)
	length*=60;
      else if (strchr(optarg,'d')!=NULL)
	length*=86400;
      break;

    case 'S':
      saltmin=atof(optarg);
      break;

    case 'A':
      altmin=atof(optarg);
      break;

    case 'd':
      dt=atof(optarg);
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

  // Error checking
  if (satno<=0) {
    fprintf(stderr,"ERROR: NORAD satellite number not provided!\n");
    return -1;
  }

  // Get TLE
  fp=fopen(tlefile,"rb");
  if (fp==NULL) {
    fprintf(stderr,"ERROR: Failed to open file with TLEs (%s)!\n",tlefile);
    return -1;
  }

  // Read TLE
  while (read_twoline(fp,satno,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);

    if (imode==SGDP4_ERROR)
      continue;

  }
  fclose(fp);

  // Object found?
  if (orb.satno!=satno) {
    fprintf(stderr,"ERROR: Object %d not found in %s!\n",satno,tlefile);
    return -1;
  }

  // Object found?
  if (orb.rev<4.0) {
    fprintf(stderr,"ERROR: Object %d is not in a LEO orbit.\n",satno);
    return -1;
  }

  // Compute radius
  if (radius<0.0)
    radius=semimajoraxis(orb);
  
  // Print header
  printf("Observer: %s (%04d) [%+.4f, %+.4f, %.0fm]\n",m.observer,m.site_id,m.lat,m.lng,m.alt*1000.0);
  printf("Elements: %s\n",tlefile);
  printf("Object: %d\n",satno);
  printf("Radius: %g km\n",radius);
  printf("Start UT Date/Time: %s for %g h \n\n",nfd,length/3600.0);
  printf("UT Date/Time         R.A.      Decl.    Azi.    Alt.  Range    Mag  Sun Alt. Type\n");
  printf("                                        (deg)   (deg) (km)          (deg)\n");
  printf("=====================================================================================\n");
  


  for (t=0.0;t<length;t+=dt) {
    // (Modified) Julian Date
    mjd=mjd0+t/86400.0;
    jd=mjd+2400000.5;
    
    // Get kepler
    tsince=1440.0*(jd-SGDP4_jd0);
    rv=sgdp4(tsince,1,&K);

    // Get positions
    obspos_xyz(mjd,&obspos,&obsvel);
    sunpos_xyz(mjd,&sunpos,&ra,&de);
    equatorial2horizontal(mjd,ra,de,&sazi,&salt);

    // Rough search first
    p0=0.0;
    p1=2.0*M_PI;
    for (i=0,pmin=0.0,mmin=15.0,altmax=0.0;i<MMAX;i++) {
      p=p0+(p1-p0)*(float) i/(float) (MMAX-1);
      state=properties(K,obspos,sunpos,radius,p,&ra,&de,&r,&mag);
      equatorial2horizontal(mjd,ra,de,&azi,&alt);

      if (opttype==0) {
	if (mag<mmin) {
	  pmin=p;
	  mmin=mag;
	}
      } else if (opttype==1) {
	if (alt>altmax && mag<15.0) {
	  pmin=p;
	  altmax=alt;
	}
      }
    }

    // Finer search
    p0=pmin-4.0*M_PI/(float) MMAX;
    p1=pmin+4.0*M_PI/(float) MMAX;
    for (i=0,pmin=0.0,mmin=15.0;i<MMAX;i++) {
      p=p0+(p1-p0)*(float) i/(float) (MMAX-1);
      state=properties(K,obspos,sunpos,radius,p,&ra,&de,&r,&mag);
      equatorial2horizontal(mjd,ra,de,&azi,&alt);

      if (opttype==0) {
	if (mag<mmin) {
	  pmin=p;
	  mmin=mag;
	}
      } else if (opttype==1) {
	if (alt>altmax && mag<15.0) {
	  pmin=p;
	  altmax=alt;
	}
      }
    }

    // Get properties before and after maximum
    pstate=properties(K,obspos,sunpos,radius,pmin-0.01,&ra,&de,&r,&mag);
    nstate=properties(K,obspos,sunpos,radius,pmin+0.01,&ra,&de,&r,&mag);
    state=properties(K,obspos,sunpos,radius,pmin,&ra,&de,&r,&mag);
    if (pstate<state && state==nstate)
      strcpy(type,"Egress");
    else if (pstate==state && state>nstate)
      strcpy(type,"Ingress");
    else if (opttype==0)
      strcpy(type,"Maximum");
    else if (opttype==1)
      strcpy(type,"Culmination");

    ra=modulo(ra,360.0);
    equatorial2horizontal(mjd,ra,de,&azi,&alt);
    azi=modulo(azi+180.0,360.0);
    mjd2date(mjd,nfd);
    dec2sex(ra/15.0,sra,0,2);
    dec2sex(de,sde,0,2);
    if (alt>altmin && salt<saltmin)
      printf("%s %s %s %6.2f %6.2f %7.1f %5.2f %6.2f   %s\n",nfd,sra,sde,azi,alt,r,mag,salt,type);
  }
  
  return 0;
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
  char abbrev[3],observer[64],filename[LIM],*env;

  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/sites.txt",env);
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"File with site information not found (%s)!\n",filename);
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
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
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

// Convert equatorial into horizontal coordinates
void equatorial2horizontal(double mjd,double ra,double de,double *azi,double *alt)
{
  double h;

  h=gmst(mjd)+m.lng-ra;
  
  *azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(m.lat*D2R)-tan(de*D2R)*cos(m.lat*D2R))*R2D,360.0);
  *alt=asin(sin(m.lat*D2R)*sin(de*D2R)+cos(m.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

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

int properties(kep_t K,xyz_t obspos,xyz_t sunpos,float radius,float t,double *ra,double *de,double *r,float *mag)
{
  double st,ct,sn,cn,si,ci;
  xyz_t satpos;
  double dx,dy,dz,rsun,rearth,psun,pearth,p;
  float phase;
  int state;
  
  // Angles
  sn=sin(K.ascn);
  cn=cos(K.ascn);
  si=sin(K.eqinc);
  ci=cos(K.eqinc);
  st=sin(K.theta+t);
  ct=cos(K.theta+t);
  
  satpos.x=-sn*ci*st+cn*ct;
  satpos.y=cn*ci*st+sn*ct;
  satpos.z=si*st;
  satpos.x*=(radius+XKMPER);
  satpos.y*=(radius+XKMPER);
  satpos.z*=(radius+XKMPER);
  
  // Sun position from satellite
  dx=-satpos.x+sunpos.x;  
  dy=-satpos.y+sunpos.y;
  dz=-satpos.z+sunpos.z;
  
  // Distances
  rsun=sqrt(dx*dx+dy*dy+dz*dz);
  rearth=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
  
  // Angles
  psun=asin(696.0e3/rsun)*R2D;
  pearth=asin(6378.135/rearth)*R2D;
  
  pearth=asin(6378.135/rearth)*R2D;
  p=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;
  p-=pearth;
  
  // Position differences
  dx=satpos.x-obspos.x;  
  dy=satpos.y-obspos.y;
  dz=satpos.z-obspos.z;
  
  // Celestial position
  *r=sqrt(dx*dx+dy*dy+dz*dz);
  *ra=atan2(dy,dx)*R2D;
  *de=asin(dz/ *r)*R2D;

  // Visibility
  if (p<-psun) {
    //    strcpy(state,"eclipsed");
    state=0;
  } else if (p>-psun && p<psun) {
    //    strcpy(state,"umbra");
    state=1;
  } else if (p>psun) {
    //    strcpy(state,"sunlit");
    state=2;
  }
  
  // Phase
  phase=acos(((obspos.x-satpos.x)*(sunpos.x-satpos.x)+(obspos.y-satpos.y)*(sunpos.y-satpos.y)+(obspos.z-satpos.z)*(sunpos.z-satpos.z))/(rsun* *r))*R2D;
  
  // Magnitude
  if (state==2) 
    *mag=STDMAG-15.0+5*log10( *r)-2.5*log10(sin(phase*D2R)+(M_PI-phase*D2R)*cos(phase*D2R));
  else
    *mag=15;

  return state;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int f,int len)
{
  int i;
  double sec,deg,min;
  char sign;
  char *form[]={"::",",,","hms","  "};

  sign=(x<0 ? '-' : ' ');
  x=3600.*fabs(x);

  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;

  if (len==7) sprintf(s,"%c%02i%c%02i%c%07.4f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==6) sprintf(s,"%c%02i%c%02i%c%06.3f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==5) sprintf(s,"%c%02i%c%02i%c%05.2f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==4) sprintf(s,"%c%02i%c%02i%c%04.1f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==2) sprintf(s,"%c%02i%c%02i%c%02i%c",sign,(int) deg,form[f][0],(int) min,form[f][1],(int) floor(sec),form[f][2]);

  return;
}
