#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qfits.h"

#define LIM 128
#define D2R M_PI/180.0
#define R2D 180.0/M_PI

struct image {
  char filename[64];
  int naxis1,naxis2;
  float *zavg;
  double ra0,de0;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  double mjd;
  float exptime;
  char nfd[32];
  int cospar;
};
struct map {
  char datadir[LIM],observer[32];
  double lng,lat;
  float alt;
  int site_id;
} m;
struct image read_fits(char *filename);

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

// Convert equatorial into horizontal coordinates
void equatorial2horizontal(double mjd,double ra,double de,double *azi,double *alt)
{
  double h;

  h=gmst(mjd)+m.lng-ra;
  
  *azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(m.lat*D2R)-tan(de*D2R)*cos(m.lat*D2R))*R2D,360.0);
  *alt=asin(sin(m.lat*D2R)*sin(de*D2R)+cos(m.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

  return;
}

int main(int argc,char *argv[])
{
  int i;
  struct image img;
  float zavg,zstd,sx,sy,wx,wy;
  double ra,de,alt,azi,mjd0=51544.5;
  char *env;
  
  // Get environment variables
  env=getenv("ST_DATADIR");
  if (env!=NULL) {
    strcpy(m.datadir,env);
  } else {
    printf("ST_DATADIR environment variable not found.\n");
  }

  // Read image
  img=read_fits(argv[1]);

  // Get site
  get_site(img.cospar);
  
  // Compute statistics
  for (i=0,zavg=0.0;i<img.naxis1*img.naxis2;i++)
    zavg+=img.zavg[i];
  zavg/=(float) img.naxis1*img.naxis2;
  for (i=0,zstd=0.0;i<img.naxis1*img.naxis2;i++)
    zstd+=pow(img.zavg[i]-zavg,2);
  zstd=sqrt(zstd/(float) (img.naxis1*img.naxis2));

  // Image scale
  sx=sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1]);
  sy=sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2]);
  wx=img.naxis1*sx/3600.0;
  wy=img.naxis2*sy/3600.0;

  // Precess to epoch of date
  precess(mjd0,img.ra0,img.de0,img.mjd,&ra,&de);

  // Get horizontal coordinates
  equatorial2horizontal(img.mjd,ra,de,&azi,&alt);
  azi=modulo(azi+180,360);
  printf("%s %14.8lf %10.6f %10.6f %10.6f %10.6f %.2f %.2f %.1f %.1f\n",argv[1],img.mjd,img.ra0,img.de0,azi,alt,img.xrms,img.yrms,zavg,zstd);

  return 0;
}

// Read fits image
struct image read_fits(char *filename)
{
  int i,j,k,l,m;
  qfitsloader ql;
  char key[FITS_LINESZ+1];
  char val[FITS_LINESZ+1];
  struct image img;
  int naxis;

  // Copy filename
  strcpy(img.filename,filename);

  // Image size
  naxis=atoi(qfits_query_hdr(filename,"NAXIS"));
  img.naxis1=atoi(qfits_query_hdr(filename,"NAXIS1"));
  img.naxis2=atoi(qfits_query_hdr(filename,"NAXIS2"));

  // MJD
  img.mjd=(double) atof(qfits_query_hdr(filename,"MJD-OBS"));
  strcpy(img.nfd,qfits_query_hdr(filename,"DATE-OBS"));
  img.exptime=atof(qfits_query_hdr(filename,"EXPTIME"));

  // COSPAR ID
  img.cospar=atoi(qfits_query_hdr(filename,"COSPAR"));

  // Transformation
  img.mjd=atof(qfits_query_hdr(filename,"MJD-OBS"));
  img.ra0=atof(qfits_query_hdr(filename,"CRVAL1"));
  img.de0=atof(qfits_query_hdr(filename,"CRVAL2"));
  img.x0=atof(qfits_query_hdr(filename,"CRPIX1"));
  img.y0=atof(qfits_query_hdr(filename,"CRPIX2"));
  img.a[0]=0.0;
  img.a[1]=3600.0*atof(qfits_query_hdr(filename,"CD1_1"));
  img.a[2]=3600.0*atof(qfits_query_hdr(filename,"CD1_2"));
  img.b[0]=0.0;
  img.b[1]=3600.0*atof(qfits_query_hdr(filename,"CD2_1"));
  img.b[2]=3600.0*atof(qfits_query_hdr(filename,"CD2_2"));
  img.xrms=3600.0*atof(qfits_query_hdr(filename,"CRRES1"));
  img.yrms=3600.0*atof(qfits_query_hdr(filename,"CRRES2"));

  // Allocate image memory
  img.zavg=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);

  // Set parameters
  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;

  // Loop over planes
  if (naxis==3)
    ql.pnum=2;
  else
    ql.pnum=0;
      
  // Initialize load
  if (qfitsloader_init(&ql) != 0) 
    printf("Error initializing data loading\n");
      
  // Test load
  if (qfits_loadpix(&ql) != 0) 
    printf("Error loading actual data\n");
  
  // Fill z array
  for (i=0,l=0;i<img.naxis1;i++) {
    for (j=0;j<img.naxis2;j++) {
      img.zavg[l]=ql.fbuf[l];
      l++;
    }
  }

  return img;
}

