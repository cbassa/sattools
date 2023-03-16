#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <wcslib/cel.h>
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

struct point {
  int flag,satno;
  double mjd,ra,de;
  float st,sr;
  char iod_line[LIM];
  xyz_t obspos;
};
struct site {
  int id;
  double lng,lat;
  float alt;
  char observer[64];
};
struct data {
  int n;
  struct point *p;
} ;
struct point decode_iod_observation(char *iod_line);
struct site get_site(int site_id);
int fgetline(FILE *file,char *s,int lim);
double modulo(double x,double y);
double gmst(double mjd);
double dgmst(double mjd);
double date2mjd(int year,int month,double day);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);
void usage();
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel);
struct data read_data(char *filename);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);

void compute_residual(char *filename,struct point p,int satno)
{
  int i,imode;
  FILE *file;
  orbit_t orb;
  xyz_t satpos,satvel;
  double dx,dy,dz;
  double r[2],ra,de;
  double rx[2],ry[2],dr,dt,drx,dry;
  double jd;
  double age;
  
  // Open catalog
  file=fopen(filename,"r");
  if (file==NULL) 
    fatal_error("Failed to open %s\n",filename);

  // Read TLE
  if (satno==0)
    read_twoline(file,p.satno,&orb);
  else
    read_twoline(file,satno,&orb);
  fclose(file);

  // Check for match
  if (orb.satno!=p.satno && satno==0) {
    //    fprintf(stderr,"object %d not found in %s\n",p.satno,filename);
    return;
  }
    
  // Initialize
  imode=init_sgdp4(&orb);
  if (imode==SGDP4_ERROR) {
    fprintf(stderr,"Error initializing SGDP4\n");
    exit(0);
  }
  
  for (i=0;i<2;i++) {
    jd=p.mjd+2400000.5+(double) i/86400;

    // Compute position
    satpos_xyz(jd,&satpos,&satvel);
    age=jd-SGDP4_jd0;

    // compute difference vector
    dx=satpos.x-p.obspos.x;  
    dy=satpos.y-p.obspos.y;
    dz=satpos.z-p.obspos.z;
  
    // Celestial position
    r[i]=sqrt(dx*dx+dy*dy+dz*dz);
    ra=modulo(atan2(dy,dx)*R2D,360.0);
    de=asin(dz/r[i])*R2D;

    // Compute offset
    forward(p.ra,p.de,ra,de,&rx[i],&ry[i]);
  }
  drx=rx[1]-rx[0];
  dry=ry[1]-ry[0];
  dt=(rx[0]*drx+ry[0]*dry)/(drx*drx+dry*dry);
  dr=(dry*rx[0]-drx*ry[0])/sqrt(drx*drx+dry*dry);
  printf("%s | %8.5f deg %9.4f sec %7.3f day, %.1f km\n",p.iod_line,dr,dt,age,r[0]);

  return;
}

void split_file(struct data d,float dtmax)
{
  int i,j,flag=0;
  FILE *file;
  char filename[LIM];
  double mjd0,dt;

  for (i=0,j=0;i<d.n;i++) {
    if (flag==1) {
      dt=86400*(d.p[i].mjd-mjd0);
      if (dt>dtmax) {
	if (file!=NULL)
	  fclose(file);
	flag=0;
	j++;
      }
    }    
    if (flag==0) {
      mjd0=d.p[i].mjd;
      flag=1;
      sprintf(filename,"split%04d.dat",j+1);
      file=fopen(filename,"w");
    }
    fprintf(file,"%s\n",d.p[i].iod_line);
  }
  if (file!=NULL)
    fclose(file);

  return;
}

int main(int argc,char *argv[])
{
  int i,arg=0,split=0;
  struct data d;
  char *datafile,catalog[LIM];
  char *env;
  float dt;
  int verbose=0,satno=0;

  env=getenv("ST_TLEDIR");
  sprintf(catalog,"%s/classfd.tle",env);
  // Decode options
  while ((arg=getopt(argc,argv,"d:c:hs:vi:"))!=-1) {
    switch(arg) {
    case 'd':
      datafile=optarg;
      break;

    case 'v':
      verbose=1;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'c':
      strcpy(catalog,optarg);
      break;

    case 's':
      dt=atof(optarg);
      split=1;
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

  // Read data
  d=read_data(datafile);

  if (verbose==1) {
    for (i=0;i<d.n;i++) {
      printf("%14.8lf %10.4f %10.4f %10.4f %10.6f %10.6f\n",d.p[i].mjd,d.p[i].obspos.x,d.p[i].obspos.y,d.p[i].obspos.z,d.p[i].ra,d.p[i].de);
    }
  }

  if (split==1) {
    split_file(d,dt);
  } else {
    for (i=0;i<d.n;i++)
      compute_residual(catalog,d.p[i],satno);
  }

  return 0;
}

// Decode IOD Observations
struct point decode_iod_observation(char *iod_line)
{
  int year,month,iday,hour,min;
  int format,epoch,me,xe,sign;
  int site_id;
  double sec,ra,mm,ss,de,dd,ds,day,mjd0;
  char secbuf[6],sn[2],degbuf[3];
  struct point p;
  struct site s;
  xyz_t vel;

  // Strip newline
  iod_line[strlen(iod_line)-1]='\0';

  // Copy full line
  strcpy(p.iod_line,iod_line);

  // Set flag
  p.flag=1;

  // Get SSN
  sscanf(iod_line,"%5d",&p.satno);

  // Get site
  sscanf(iod_line+16,"%4d",&site_id);
  s=get_site(site_id);

  // Skip if site not found
  if (s.id<0) {
    fprintf(stderr,"Site %d not found!\n",site_id);
    p.flag=0;
  }
  
  // Decode date/time
  sscanf(iod_line+23,"%4d%2d%2d%2d%2d%5s",&year,&month,&iday,&hour,&min,secbuf);
  sec=atof(secbuf);
  sec/=pow(10,strlen(secbuf)-2);
  day=(double) iday+(double) hour/24.0+(double) min/1440.0+(double) sec/86400.0;
  p.mjd=date2mjd(year,month,day);

  // Get uncertainty in time
  sscanf(iod_line+41,"%1d%1d",&me,&xe);
  p.st=(float) me*pow(10,xe-8);

  // Get observer position
  obspos_xyz(p.mjd,s.lng,s.lat,s.alt,&p.obspos,&vel);

  // Skip empty observations
  if (strlen(iod_line)<64 || (iod_line[54]!='+' && iod_line[54]!='-')) 
    p.flag=0;

  // Get format, epoch
  sscanf(iod_line+44,"%1d%1d",&format,&epoch);

  // Read position
  sscanf(iod_line+47,"%2lf%2lf%3lf%1s",&ra,&mm,&ss,sn);
  sscanf(iod_line+55,"%2lf%2lf%2s",&de,&dd,degbuf);
  ds=atof(degbuf);
  if (strlen(degbuf)==1)
    ds*=10;
  sign=(sn[0]=='-') ? -1 : 1;
  sscanf(iod_line+62,"%1d%1d",&me,&xe);
  p.sr=(float) me*pow(10,xe-8);
  
  // Decode position
  switch(format) 
    {
      // Format 1: RA/DEC = HHMMSSs+DDMMSS MX   (MX in seconds of arc)
    case 1 : 
      ra+=mm/60+ss/36000;
      de=sign*(de+dd/60+ds/3600);
      p.sr/=3600.0;
      break;
      // Format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
    case 2:
      ra+=mm/60+ss/60000;
      de=sign*(de+dd/60+ds/6000);
      p.sr/=60.0;
      break;
      // Format 3: RA/DEC = HHMMmmm+DDdddd MX   (MX in degrees of arc)
    case 3 :
      ra+=mm/60+ss/60000;
      de=sign*(de+dd/100+ds/10000);
      break;
      // Format 7: RA/DEC = HHMMSSs+DDdddd MX   (MX in degrees of arc)
    case 7 :
      ra+=mm/60+ss/36000;
      de=sign*(de+dd/100+ds/10000);
      break;
    default :
      fprintf(stderr,"IOD Format not implemented\n");
      p.flag=0;
      break;
  }
  // Convert to degrees
  ra*=15.0;

  // Get precession epoch
  if (epoch==0) {
    p.ra=ra;
    p.de=de;
    return p;
  } else if (epoch==4) {
    mjd0=33281.9235;
  } else if (epoch==5) {
    mjd0=51544.5;
  } else {
    fprintf(stderr,"Observing epoch not implemented\n");
    p.flag=0;
  }

  // Precess position
  precess(mjd0,ra,de,p.mjd,&p.ra,&p.de);

  return p;
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
    return s;
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

  if (id!=site_id)
    s.id==-1;

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

void usage()
{
  printf("bla\n");

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

// Read data
struct data read_data(char *filename)
{
  int i=0;
  char line[LIM];
  FILE *file;
  struct data d;

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
    if (isdigit(line[0]))
      d.p[i++]=decode_iod_observation(line);
  }

  // Close file
  fclose(file);

  return d;
}

// Get a x and y from a RA and Decl
void forward(double ra0,double de0,double ra,double de,double *x,double *y)
{
  int i,status;
  double phi,theta;
  struct celprm cel;

  // Initialize Reference Angles
  celini(&cel);
  cel.ref[0]=ra0;
  cel.ref[1]=de0;
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;
  strcpy(cel.prj.code,"TAN");

  if (celset(&cel)) {
    printf("Error in Projection (celset)\n");
    return;
  }
  cels2x(&cel,1,0,1,1,&ra,&de,&phi,&theta,x,y,&status);

  return;
}
