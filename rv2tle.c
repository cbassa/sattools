#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
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

// Dot product
float dot(xyz_t a,xyz_t b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}

// Magnitude
double magnitude(xyz_t a)
{
  return sqrt(dot(a,a));
}

// Cross product
xyz_t cross(xyz_t a,xyz_t b)
{
  xyz_t c;

  c.x=a.y*b.z-a.z*b.y;
  c.y=a.z*b.x-a.x*b.z;
  c.z=a.x*b.y-a.y*b.x;

  return c;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
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

// Clasical elements
orbit_t classel(int ep_year,double ep_day,xyz_t r,xyz_t v)
{
  int i;
  double rm,vm,vm2,rvm,mu=1.0;;
  double chi,xp,yp,sx,cx,b,ee;
  double a,ecc,incl,node,peri,mm,n;
  xyz_t h,e,kk,nn;
  orbit_t orb;

  r.x/=XKMPER;
  r.y/=XKMPER;
  r.z/=XKMPER;
  v.x/=(XKE*XKMPER/AE*XMNPDA/86400.0);
  v.y/=(XKE*XKMPER/AE*XMNPDA/86400.0);
  v.z/=(XKE*XKMPER/AE*XMNPDA/86400.0);

  rm=magnitude(r);
  vm2=dot(v,v);
  rvm=dot(r,v);
  h=cross(r,v);
  chi=dot(h,h)/mu;

  e.x=(vm2/mu-1.0/rm)*r.x-rvm/mu*v.x;
  e.y=(vm2/mu-1.0/rm)*r.y-rvm/mu*v.y;
  e.z=(vm2/mu-1.0/rm)*r.z-rvm/mu*v.z;

  a=pow(2.0/rm-vm2/mu,-1);
  ecc=magnitude(e);
  incl=acos(h.z/magnitude(h))*R2D;
  
  kk.x=0.0;
  kk.y=0.0;
  kk.z=1.0;
  nn=cross(kk,h);
  if (nn.x==0.0 && nn.y==0.0)
    nn.x=1.0;
  node=atan2(nn.y,nn.x)*R2D;
  if (node<0.0)
    node+=360.0;

  peri=acos(dot(nn,e)/(magnitude(nn)*ecc))*R2D;
  if (e.z<0.0)
    peri=360.0-peri;
  if (peri<0.0)
    peri+=360.0;

  // Elliptic motion
  if (ecc<1.0) {
    xp=(chi-rm)/ecc;
    yp=rvm/ecc*sqrt(chi/mu);
    b=a*sqrt(1.0-ecc*ecc);
    cx=xp/a+ecc;
    sx=yp/b;
    ee=atan2(sx,cx);
    n=XKE*sqrt(mu/(a*a*a));
    mm=(ee-ecc*sx)*R2D;
  } 
  if (mm<0.0)
    mm+=360.0;

  // Fill
  orb.satno=0;
  orb.eqinc=incl*D2R;
  orb.ascn=node*D2R;
  orb.argp=peri*D2R;
  orb.mnan=mm*D2R;
  orb.ecc=ecc;
  orb.rev=XKE*pow(a,-3.0/2.0)*XMNPDA/(2.0*M_PI);
  orb.bstar=0.0;
  orb.ep_year=ep_year;
  orb.ep_day=ep_day;
  orb.norb=0;

  return orb;
}

orbit_t rv2el(int satno,double mjd,xyz_t r0,xyz_t v0)
{
  int i,imode;
  orbit_t orb[5],orb1[5];
  xyz_t r,v;
  kep_t kep;
  char line1[70],line2[70];
  int ep_year;
  double ep_day;

  // Epoch
  ep_day=mjd2doy(mjd,&ep_year);

  // Initial guess
  orb[0]=classel(ep_year,ep_day,r0,v0);
  orb[0].satno=satno;
  
  for (i=0;i<5;i++) {
    // Propagate
    imode=init_sgdp4(&orb[i]);
    imode=satpos_xyz(mjd+2400000.5,&r,&v);

    // Compute initial orbital elements
    orb1[i]=classel(ep_year,ep_day,r,v);

    // Adjust
    orb[i+1].rev=orb[i].rev+orb[0].rev-orb1[i].rev;
    orb[i+1].ascn=orb[i].ascn+orb[0].ascn-orb1[i].ascn;
    orb[i+1].argp=orb[i].argp+orb[0].argp-orb1[i].argp;
    orb[i+1].mnan=orb[i].mnan+orb[0].mnan-orb1[i].mnan;
    orb[i+1].eqinc=orb[i].eqinc+orb[0].eqinc-orb1[i].eqinc;
    orb[i+1].ecc=orb[i].ecc+orb[0].ecc-orb1[i].ecc;
    orb[i+1].ep_year=orb[i].ep_year;
    orb[i+1].ep_day=orb[i].ep_day;
    orb[i+1].satno=orb[i].satno;
    orb[i+1].norb=orb[i].norb;
    orb[i+1].bstar=orb[i].bstar;

    // Keep in range
    if (orb[i+1].ecc<0.0)
      orb[i+1].ecc=0.0;
    if (orb[i+1].eqinc<0.0)
      orb[i+1].eqinc=0.0;
    orb[i+1].mnan=modulo(orb[i+1].mnan,2.0*M_PI);
    orb[i+1].ascn=modulo(orb[i+1].ascn,2.0*M_PI);
    orb[i+1].argp=modulo(orb[i+1].argp,2.0*M_PI);
  }

  return orb[i];
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

void usage(void)
{

  return;
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

int main(int argc,char *argv[])
{
  int imode,satno=99999,arg;
  FILE *file;
  orbit_t orb;
  xyz_t r,v;
  char xyzfile[LIM],nfd[32],line[LIM];
  char line1[80],line2[80],desig[20]="14999A";
  double mjd;
  char *env;

  // Decode options
  while ((arg=getopt(argc,argv,"p:i:d:"))!=-1) {
    switch (arg) {

    case 'p':
      strcpy(xyzfile,optarg);
      break;

    case 'i':
      satno=atoi(optarg);
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
  file=fopen(xyzfile,"r");
  while (fgetline(file,line,LIM)>0) {
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf",&mjd,&r.x,&r.y,&r.z,&v.x,&v.y,&v.z);
      
    // Convert
    orb=rv2el(orb.satno,mjd,r,v);
    orb.satno=satno;
    strcpy(orb.desig,desig);

    format_tle(orb,line1,line2);
    printf("%s\n%s\n",line1,line2);
  }
  fclose(file);

  return 0;
}

