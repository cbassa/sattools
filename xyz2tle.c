#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "sgdp4h.h"
#include "satutl.h"


#define LIM 80
#define NMAX 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKE 0.07436680 // Guassian Gravitational Constant
#define XKMPER 6378.135
#define AE 1.0
#define XMNPDA 1440.0

struct data {
  int n,nsel;
  struct point *p;
  double chisq,rms;
} d;
struct point {
  int flag;
  double mjd;
  xyz_t r;
};
orbit_t orb;

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

// Read data file
struct data read_data(char *filename,double mjd0)
{
  int i=0,status;
  char line[LIM];
  FILE *file;
  struct data d;
  int min;
  double ra,de,ra0,de0,r;
  double x,y,z;

  // Open file
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"Failed to open %s\n",filename);
    exit(1);
  }

  // Count lines
  while (fgetline(file,line,LIM)>0) 
    i++;
  d.n=i;

  // Allocate
  d.p=(struct point *) malloc(sizeof(struct point)*d.n);

  // Rewind file
  rewind(file);

  // Read data
  i=0;
  while (fgetline(file,line,LIM)>0) {
    status=sscanf(line,"%d,%lf,%lf,%lf",&min,&x,&y,&z);
    d.p[i].mjd=mjd0+(double) min/1440.0;

    // Precess position
    r=sqrt(x*x+y*y+z*z);
    ra0=atan2(y,x)*R2D;
    de0=asin(z/r)*R2D;

    precess(51544.5,ra0,de0,d.p[i].mjd,&ra,&de);
    d.p[i].r.x=r*cos(de*D2R)*cos(ra*D2R);
    d.p[i].r.y=r*cos(de*D2R)*sin(ra*D2R);
    d.p[i].r.z=r*sin(de*D2R);

    d.p[i].flag=0;
    i++;
  }

  // Close file
  fclose(file);

  return d;
}

// Read tle
orbit_t read_tle(char *filename,int satno)
{
  int i;
  FILE *file;
  orbit_t orb;

  file=fopen(filename,"r");
  if (file==NULL) 
    fatal_error("Failed to open %s\n",filename);

  // Read TLE
  read_twoline(file,satno,&orb);
  fclose(file);

  return orb;
}

// Chi-squared
double chisq(double *a)
{
  int i,imode,n;
  double chisq;
  xyz_t satpos,satvel,dr;

  // Construct struct
  // a[0]: inclination
  // a[1]: RA of ascending node
  // a[2]: eccentricity
  // a[3]: argument of periastron
  // a[4]: mean anomaly
  // a[5]: revs per day
  // a[6]: bstar drag

  if (a[2]<0.0)
    a[2]=0.0;
  if (a[0]<0.0) {
    a[0]*=-1;
    a[1]+=180.0;
  } else if (a[0]>180.0) {
    a[0]=180.0;
  }
  if (a[5]>20.0)
    a[5]=20.0;
  if (a[5]<0.1)
    a[5]=0.1;

  // Set parameters
  orb.eqinc=RAD(a[0]);
  orb.ascn=RAD(modulo(a[1],360.0));
  orb.ecc=a[2];
  orb.argp=RAD(modulo(a[3],360.0));
  orb.mnan=RAD(modulo(a[4],360.0));
  orb.rev=a[5];
  orb.bstar=a[6];

  // Initialize
  imode=init_sgdp4(&orb);
  if (imode==SGDP4_ERROR)
    printf("Error\n");

  // Loop over points
  for (i=0,chisq=0.0,n=0;i<d.n;i++) {
    // Skip unflagged positions
    if (d.p[i].flag!=1)
      continue;

    // Get satellite position
    satpos_xyz(d.p[i].mjd+2400000.5,&satpos,&satvel);

    dr.x=(satpos.x-d.p[i].r.x);
    dr.y=(satpos.y-d.p[i].r.y);
    dr.z=(satpos.z-d.p[i].r.z);

    // Add
    chisq+=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
    n++;
  }
  chisq/=(double) n;

  d.rms=sqrt(chisq);
  d.nsel=n;

  return chisq;
}

double decode_filename(char *filename,int *satno)
{
  int year,month,day,hour,min,sec;
  int status;
  double mjd;
  
  status=sscanf(filename,"%6d_%4d%2d%2d_%2d%2d%2d",satno,&year,&month,&day,&hour,&min,&sec);

  mjd=date2mjd(year,month,(double) day+hour/24.0+min/1440.0+sec/86400.0);

  return mjd;
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

// State vector to SGP4 elements
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
  
  for (i=0;i<4;i++) {
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
  }

  return orb[i];
}

// Fit
void fit(orbit_t orb,int *ia)
{
  int i,n;
  double a[7],da[7];
  double db[7]={0.1,0.1,0.002,0.1,0.1,0.01,0.0001};

  // Copy parameters
  a[0]=orb.eqinc*R2D;
  da[0]=da[0]*R2D;
  a[1]=orb.ascn*R2D;
  da[1]=da[1]*R2D;
  a[2]=orb.ecc;
  a[3]=orb.argp*R2D;
  da[3]=da[3]*R2D;
  a[4]=orb.mnan*R2D;
  da[4]=da[4]*R2D;
  a[5]=orb.rev;
  a[6]=orb.bstar;

  for (i=0;i<7;i++) {
    if (ia[i]==1)
      da[i]=db[i];
    else
      da[i]=0.0;
  }

  // Construct struct
  // a[0]: inclination
  // a[1]: RA of ascending node
  // a[2]: eccentricity
  // a[3]: argument of periastron
  // a[4]: mean anomaly
  // a[5]: revs per day
  // a[6]: bstar

  // Count highlighted points
  for (i=0,n=0;i<d.n;i++)
    if (d.p[i].flag==1)
      n++;

  if (n>0)
    versafit(n,7,a,da,chisq,0.0,1e-7,"n");

  // Return parameters
  orb.eqinc=RAD(a[0]);
  orb.ascn=RAD(modulo(a[1],360.0));
  orb.ecc=a[2];
  orb.argp=RAD(modulo(a[3],360.0));
  orb.mnan=RAD(modulo(a[4],360.0));
  orb.rev=a[5];
  orb.bstar=a[6];

  return;
}

int main(int argc,char *argv[])
{
  int i,j,k,arg=0,satno=0,satname=0,usecatalog=0,imode,m=10;
  char *datafile,*catalog,filename[32];
  int ia[7]={0,0,0,0,0,0,0};
  char line1[70],line2[70],desig[10];
  double mjd;
  xyz_t r,v;
  FILE *file;

  // Decode options
  while ((arg=getopt(argc,argv,"d:c:i:"))!=-1) {
    switch(arg) {

    case 'd':
      datafile=optarg;
      break;

    case 'c':
      catalog=optarg;
      usecatalog=1;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    default:
      return 0;
    }
  }

  // Reloop stderr
  freopen("/tmp/stderr.txt","w",stderr);

  // Decode filename
  mjd=decode_filename(datafile,&satname);

  // Read data
  d=read_data(datafile,mjd);

  // Write data
  sprintf(filename,"%06d.xyz",satname);
  file=fopen(filename,"w");
  for (i=0;i<d.n;i+=10) 
    fprintf(file,"%lf %f %f %f\n",d.p[i].mjd,d.p[i].r.x,d.p[i].r.y,d.p[i].r.z);
  fclose(file);

  // Open elements
  sprintf(filename,"%06d.tle",satname);
  file=fopen(filename,"w");

  // Estimate orbit
  k=5040;
  if (usecatalog==0) {
    // Set initial state vector
    r.x=d.p[k].r.x;
    r.y=d.p[k].r.y;
    r.z=d.p[k].r.z;
    v.x=(d.p[k+1].r.x-d.p[k].r.x)/60.0;
    v.y=(d.p[k+1].r.y-d.p[k].r.y)/60.0;
    v.z=(d.p[k+1].r.z-d.p[k].r.z)/60.0;
    
    // Estimate initial orbit from elements
    orb=rv2el(99999,d.p[k].mjd,r,v);
    strcpy(orb.desig,"14999A");
  } else {
    // Read orbit
    orb=read_tle(catalog,satno);
    strcpy(desig,orb.desig);
    // Propagate
    imode=init_sgdp4(&orb);
    imode=satpos_xyz(d.p[k].mjd+2400000.5,&r,&v);
    orb=rv2el(orb.satno,d.p[k].mjd,r,v);
    //      orb.satno=99999;
    //      strcpy(orb.desig,"14999A");
    strcpy(orb.desig,desig);
  }

  // Set flags
  for (j=0;j<d.n;j++)
    d.p[j].flag=0;
    
  for (j=0;j<d.n;j+=10)
    d.p[j].flag=1;
  
  // Fit orbit
  for (j=0;j<10;j++) {
    if (j==1) ia[4]=1;
    if (j==2) ia[1]=1;
    if (j==3) ia[0]=1;
    if (j==4) ia[5]=1;
    if (j==5) ia[3]=1;
    if (j==6) ia[2]=1;
    if (j==7) ia[6]=1;
    fit(orb,ia);
  }

  // Format TLE
  format_tle(orb,line1,line2);
  fprintf(file,"SO %d\n%s\n%s\n# %d positions, %.1f km rms\n",satname,line1,line2,d.nsel,d.rms);
  printf("SO %d\n%s\n%s\n# %d positions, %.1f km rms\n",satname,line1,line2,d.nsel,d.rms);

  // Close output file
  fclose(file);

  return 0;
}
