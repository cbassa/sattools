#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <wcslib/cel.h>
#include <cpgplot.h>
#include "qfits.h"
#include "sgdp4h.h"

#define LIM 80
#define MMAX 10
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)
#define STDMAG 6.0

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

struct map {
  double lat,lng;
  float alt;
  char observer[32];
  int site_id;
} m;
struct point {
  double mjd;
  xyz_t obspos,sunpos;
  double zeta,z,theta;
} p[MMAX];
struct image {
  char filename[64];
  int naxis,naxis1,naxis2,nframes;
  float *zavg,*zstd,*zmax,*znum;
  double ra0,de0;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  double mjd;
  float *dt,exptime;
  char nfd[32];
  int cospar,tracked;
};
struct image read_fits(char *filename);
void get_site(int site_id);
double modulo(double,double);
void obspos_xyz(double,xyz_t *,xyz_t *);
void sunpos_xyz(double,xyz_t *);
double gmst(double);
double dgmst(double);
void precession_angles(double mjd0,double mjd,double *zeta,double *z,double *theta);
void initialize(struct image img);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double ra0,double de0,double x,double y,double *ra,double *de);
struct sat apparent_position(double mjd);

void plot_satellites(char *tlefile,struct image img,long satno,double mjd0,float dt,int color)
{
  int i;
  orbit_t orb;
  int imode,flag,textflag,sflag;
  FILE *fp=NULL,*file;
  xyz_t satpos,satvel;
  float x,y,x0,y0;
  char norad[7],satname[30],state[16];
  float isch;
  char filename[128];
  double dx,dy,dz,r,ra,de,d,rsun,rearth;
  double psun,pearth,ptot;
  double a,b,c;
  double rx,ry;

  cpgqch(&isch);

  // Image determinant
  d=img.a[1]*img.b[2]-img.a[2]*img.b[1];

  // Open TLE file
  fp=fopen(tlefile,"rb");
  if (fp==NULL)
    return;

  cpgsci(color);

  // Open file
  sprintf(filename,"%s.id",img.filename);
  file=fopen(filename,"a");

  // Read TLEs
  while (read_twoline(fp,satno,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);

    sprintf(norad," %05ld",Isat);

    if (imode==SGDP4_ERROR)
      continue;

    // Loop over times
    for (flag=0,textflag=0,sflag=0,i=0;i<MMAX;i++) {
      // Satellite position
      satpos_xyz(p[i].mjd+2400000.5,&satpos,&satvel);

      // Check on radius
      r=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
      if (r>300000)
      	continue;
      
      // Relative to observer
      dx=satpos.x-p[i].obspos.x;  
      dy=satpos.y-p[i].obspos.y;
      dz=satpos.z-p[i].obspos.z;

      // Celestial position
      r=sqrt(dx*dx+dy*dy+dz*dz);
      ra=modulo(atan2(dy,dx),2.0*M_PI);
      de=asin(dz/r);

      
      // Correct for precession
      a=cos(de)*sin(ra+p[i].zeta);
      b=cos(p[i].theta)*cos(de)*cos(ra+p[i].zeta)-sin(p[i].theta)*sin(de);
      c=sin(p[i].theta)*cos(de)*cos(ra+p[i].zeta)+cos(p[i].theta)*sin(de);
      ra=modulo((atan2(a,b)+p[i].z)*R2D,360.0);
      de=asin(c)*R2D;
      
      // Adjust for stationary camera
      if (img.tracked==0)
	ra+=gmst(img.mjd+0.5*img.exptime/86400.0)-gmst(p[i].mjd);
      
      // Check if nearby enough
      r=acos(sin(img.de0*D2R)*sin(de*D2R)+cos(img.de0*D2R)*cos(de*D2R)*cos((img.ra0-ra)*D2R))*R2D;
      if (r<90.0)
	forward(img.ra0,img.de0,ra,de,&rx,&ry);
      else
	continue;

      // Satellite position relative to the Sun
      dx=-satpos.x+p[i].sunpos.x;  
      dy=-satpos.y+p[i].sunpos.y;
      dz=-satpos.z+p[i].sunpos.z;

      // Distances
      rsun=sqrt(dx*dx+dy*dy+dz*dz);
      rearth=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);

      // Angles
      psun=asin(696.0e3/rsun)*R2D;
      pearth=asin(6378.135/rearth)*R2D;
      ptot=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;

      // Visibility state
      if (ptot-pearth<-psun) {
	cpgsls(4);
	strcpy(state,"eclipsed");
      } else if (ptot-pearth>-psun && ptot-pearth<psun) {
	cpgsls(2);
	strcpy(state,"umbra");
	sflag=1;
      } else if (ptot-pearth>psun) {
	cpgsls(1);
	strcpy(state,"sunlit");
	sflag=1;
      }

      // Convert image position
      dx=rx-img.a[0];
      dy=ry-img.b[0];
      x=(img.b[2]*dx-img.a[2]*dy)/d+img.x0;
      y=(img.a[1]*dy-img.b[1]*dx)/d+img.y0;

      // Print name if in viewport
      if (x>0.0 && x<img.naxis1 && y>0.0 && y<img.naxis2 && textflag==0) {
	if (flag!=0)
	  cpgdraw(x,y);
	cpgsch(0.65);
	cpgtext(x,y,norad);
	cpgsch(isch);
	cpgmove(x,y);
	textflag=1;
      }
      if (i==0) {
	x0=x;
	y0=y;
      }

      // Plot satellites
      if (flag==0) {
	cpgpt1(x,y,17);
	cpgmove(x,y);
	flag=1;
      } else {
	cpgdraw(x,y);
      }
    }
    if (textflag==1) {
      if (sflag==0)
	fprintf(file,"%.23s %8.3f %8.3f %8.3f %8.3f %8.5f %s %s eclipsed\n",img.nfd+1,x0,y0,x,y,img.exptime,norad,tlefile);
      else
	fprintf(file,"%.23s %8.3f %8.3f %8.3f %8.3f %8.5f %s %s sunlit\n",img.nfd+1,x0,y0,x,y,img.exptime,norad,tlefile);
    }
  }
  fclose(fp);
  fclose(file);
  cpgsci(1); 

  return;
}


int main(int argc,char *argv[])
{
  int i;
  struct image img;
  float zmin,zmax,zavg,zstd;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  char text[128];
  char *env,filename[128];
  float sx,sy,wx,wy;
  
  // Read fits file
  img=read_fits(argv[1]);

  // Set site
  get_site(img.cospar);

  // Initialize
  initialize(img);

  // Fill buffer
  if (img.naxis==3) {
    for (i=0,zavg=0.0;i<img.naxis1*img.naxis2;i++)
      zavg+=img.zmax[i];
    zavg/=(float) img.naxis1*img.naxis2;
    for (i=0,zstd=0.0;i<img.naxis1*img.naxis2;i++)
      zstd+=pow(img.zmax[i]-zavg,2);
    zstd=sqrt(zstd/(float) (img.naxis1*img.naxis2));
    zmin=zavg-2*zstd;
    zmax=zavg+6*zstd;
  } else {
    for (i=0,zavg=0.0;i<img.naxis1*img.naxis2;i++)
      zavg+=img.zavg[i];
    zavg/=(float) img.naxis1*img.naxis2;
    for (i=0,zstd=0.0;i<img.naxis1*img.naxis2;i++)
      zstd+=pow(img.zavg[i]-zavg,2);
    zstd=sqrt(zstd/(float) (img.naxis1*img.naxis2));
    zmin=zavg-2*zstd;
    zmax=zavg+6*zstd;
  }

  // Sizes
  sx=sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1]);
  sy=sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2]);
  wx=img.naxis1*sx/3600.0;
  wy=img.naxis2*sy/3600.0;
  
  if (argc==3)
    cpgopen(argv[2]);
  else
    cpgopen("/xs");
  cpgpap(14.,1.0);
  cpgsvp(0.1,0.95,0.1,0.8);

  cpgsch(0.8);
  sprintf(text,"UT Date: %.23s  COSPAR ID: %04d",img.nfd+1,img.cospar);
  cpgmtxt("T",6.0,0.0,0.0,text);
  sprintf(text,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')",img.ra0,img.xrms,img.de0,img.yrms);
  if (img.xrms<1e-3 || img.yrms<1e-3 || img.xrms/sx>2.0 || img.yrms/sy>2.0)
    cpgsci(2);
  else
    cpgsci(1);
  cpgmtxt("T",4.8,0.0,0.0,text);
  cpgsci(1);

  
  sprintf(text,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d",wx,wy,sx,sy);
  cpgmtxt("T",3.6,0.0,0.0,text);
  sprintf(text,"Stat: %5.1f+-%.1f (%.1f-%.1f)",zavg,zstd,zmin,zmax);
  cpgmtxt("T",2.4,0.0,0.0,text);
  
  cpgsch(1.0);
  cpgwnad(0.0,img.naxis1,0.0,img.naxis2);

  cpglab("x (pix)","y (pix)"," ");
  cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
    
  if (img.naxis==3)
    cpgimag(img.zmax,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
  else
    cpgimag(img.zavg,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
  cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

  cpgstbg(1);

  // Environment variables
  env=getenv("ST_TLEDIR");
  sprintf(filename,"%s/classfd.tle",env);
  plot_satellites(filename,img,0,img.mjd,img.exptime,4);
  sprintf(filename,"%s/inttles.tle",env);
  plot_satellites(filename,img,0,img.mjd,img.exptime,3);
  sprintf(filename,"%s/catalog.tle",env);
  plot_satellites(filename,img,0,img.mjd,img.exptime,0);
  sprintf(filename,"%s/jsc.txt",env);
  plot_satellites(filename,img,0,img.mjd,img.exptime,5);

  cpgend();

  
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

  // Copy filename
  strcpy(img.filename,filename);

  // Image size
  img.naxis=atoi(qfits_query_hdr(filename,"NAXIS"));
  img.naxis1=atoi(qfits_query_hdr(filename,"NAXIS1"));
  img.naxis2=atoi(qfits_query_hdr(filename,"NAXIS2"));

  // MJD
  img.mjd=(double) atof(qfits_query_hdr(filename,"MJD-OBS"));
  strcpy(img.nfd,qfits_query_hdr(filename,"DATE-OBS"));
  img.exptime=atof(qfits_query_hdr(filename,"EXPTIME"));

  // COSPAR ID
  img.cospar=atoi(qfits_query_hdr(filename,"COSPAR"));

  // Tracked
  if (qfits_query_hdr(filename,"TRACKED")!=NULL)
    img.tracked=atoi(qfits_query_hdr(filename,"TRACKED"));
  else
    img.tracked=0;

  // Transformation
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

  // Set parameters
  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;

  // Read four-frame info
  if (img.naxis==3) {
    // Number of frames
    img.nframes=atoi(qfits_query_hdr(filename,"NFRAMES"));

    // Timestamps
    img.dt=(float *) malloc(sizeof(float)*img.nframes);
    for (i=0;i<img.nframes;i++) {
      sprintf(key,"DT%04d",i);
      strcpy(val,qfits_query_hdr(filename,key));
      sscanf(val+1,"%f",&img.dt[i]);
    //    img.dt[i]=atof(qfits_query_hdr(filename,key));
    }

    // Allocate image memory
    img.zavg=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
    img.zstd=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
    img.zmax=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
    img.znum=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);

    // Loop over planes
    for (k=0;k<4;k++) {
      ql.pnum=k;
      
      // Initialize load
      if (qfitsloader_init(&ql) != 0) 
	printf("Error initializing data loading\n");
      
      // Test load
      if (qfits_loadpix(&ql) != 0) 
	printf("Error loading actual data\n");
      
      // Fill z array
      for (i=0,l=0;i<img.naxis1;i++) {
	for (j=0;j<img.naxis2;j++) {
	  if (k==0) img.zavg[l]=ql.fbuf[l];
	  if (k==1) img.zstd[l]=ql.fbuf[l];
	  if (k==2) img.zmax[l]=ql.fbuf[l];
	  if (k==3) img.znum[l]=ql.fbuf[l];
	  l++;
	}
      }
    }
  } else {
    // Allocate image memory
    img.zavg=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
    
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
  }

  return img;
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
void sunpos_xyz(double mjd,xyz_t *pos)
{
  double jd,t,l0,m,e,c,r;
  double n,s,ecl,ra,de;

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

  ra=atan2(cos(ecl)*sin(s),cos(s));
  de=asin(sin(ecl)*sin(s));

  pos->x=r*cos(de)*cos(ra)*XKMPAU;
  pos->y=r*cos(de)*sin(ra)*XKMPAU;
  pos->z=r*sin(de)*XKMPAU;

  return;
}

// Compute precession angles
void precession_angles(double mjd0,double mjd,double *zeta,double *z,double *theta)
{
  double t0,t;

  // Time in centuries
  t0=(mjd0-51544.5)/36525.0;
  t=(mjd-mjd0)/36525.0;

  // Precession angles
  *zeta=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *zeta+=(0.30188-0.000344*t0)*t*t+0.017998*t*t*t;
  *zeta*=D2R/3600.0;
  *z=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *z+=(1.09468+0.000066*t0)*t*t+0.018203*t*t*t;
  *z*=D2R/3600.0;
  *theta=(2004.3109-0.85330*t0-0.000217*t0*t0)*t;
  *theta+=-(0.42665+0.000217*t0)*t*t-0.041833*t*t*t;
  *theta*=D2R/3600.0;
  
  return;
}

// Initialize observer and sun position and precession angles
void initialize(struct image img)
{
  int i;
  double t;
  xyz_t obsvel;

  // Loop over points
  for (i=0;i<MMAX;i++) {
    // Compute time
    t=img.exptime*(float) i/(float) (MMAX-1);
    p[i].mjd=img.mjd+t/86400.0;

    // Compute observer and sun position
    obspos_xyz(p[i].mjd,&p[i].obspos,&obsvel);
    sunpos_xyz(p[i].mjd,&p[i].sunpos);

    // Compute precession angles
    precession_angles(p[i].mjd,51544.5,&p[i].zeta,&p[i].z,&p[i].theta);
  }
    
  return;
}

