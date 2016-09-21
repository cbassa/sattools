#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cel.h"
#include "cpgplot.h"
#include "qfits.h"
#include "sgdp4h.h"

#define LIM 80
#define NMAX 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)
#define STDMAG 6.0
#define MMAX 10

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

struct map {
  double lat,lng;
  float alt;
  char observer[32];
  int site_id;
} m;
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
struct sat {
  long Isat;
  char state[10];
  float mag;
  double jd;
  double dx,dy,dz;
  double x,y,z,vx,vy,vz;
  double rsun,rearth;
  double psun,pearth,p,phase;
  double r,v,ra,de;
  double azi,alt;
  double rx,ry;
};
struct image read_fits(char *filename);
struct sat apparent_position(double mjd);
double modulo(double,double);
void obspos_xyz(double,xyz_t *,xyz_t *);
void sunpos_xyz(double,xyz_t *);
double gmst(double);
double dgmst(double);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double ra0,double de0,double x,double y,double *ra,double *de);

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

void plot_satellites(char *tlefile,struct image img,long satno,double mjd0,float dt,int color)
{
  int i;
  orbit_t orb;
  struct sat s;
  int imode,flag,textflag;
  FILE *fp=NULL,*file;;
  xyz_t satpos,obspos,satvel,sunpos;
  double mjd,jd,dx,dy,dz;
  double rx,ry,ra,de,azi,alt,r,t,d;
  float x,y,x0,y0;
  char norad[7],satname[30];
  float isch;
  float rsun,rearth,psun,pearth,p;
  char filename[128];

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

    for (flag=0,textflag=0,i=0;i<MMAX;i++) {
      t=img.exptime*(float) i/(float) (MMAX-1);
      mjd=mjd0+t/86400.0;

      // Compute apparent position
      s=apparent_position(mjd);

      // Adjust for stationary camera
      if (img.tracked==0) 
	s.ra+=gmst(img.mjd+0.5*img.exptime/86400.0)-gmst(mjd);

      // Convert to rx,ry
      r=acos(sin(img.de0*D2R)*sin(s.de*D2R)+cos(img.de0*D2R)*cos(s.de*D2R)*cos((img.ra0-s.ra)*D2R))*R2D;
      if (r<90.0)
	forward(img.ra0,img.de0,s.ra,s.de,&s.rx,&s.ry);
      else
	continue;

      // Convert image position
      dx=s.rx-img.a[0];
      dy=s.ry-img.b[0];
      x=(img.b[2]*dx-img.a[2]*dy)/d+img.x0;
      y=(img.a[1]*dy-img.b[1]*dx)/d+img.y0;

      // Visibility
      if (s.p-s.pearth<-s.psun) {
	cpgsls(4);
      } else if (s.p-s.pearth>-s.psun && s.p-s.pearth<s.psun) {
	cpgsls(2);
      } else if (s.p-s.pearth>s.psun) {
	cpgsls(1);
      }

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
    if (textflag==1)
      fprintf(file,"%.23s %8.3f %8.3f %8.3f %8.3f %8.5f %s %s\n",img.nfd+1,x0,y0,x,y,img.exptime,norad,tlefile);

  }
  fclose(fp);
  fclose(file);
  cpgsci(1); 

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

  img=read_fits(argv[1]);

  // Set site
  get_site(img.cospar);

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

  if (argc==3)
    cpgopen(argv[2]);
  else
    cpgopen("/xs");
  cpgpap(0.,1.0);
  cpgsvp(0.1,0.95,0.1,0.8);

  cpgsch(0.8);
  sprintf(text,"UT Date: %.23s  COSPAR ID: %04d",img.nfd+1,img.cospar);
  cpgmtxt("T",6.0,0.0,0.0,text);
  sprintf(text,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')",img.ra0,img.xrms,img.de0,img.yrms);
  cpgmtxt("T",4.8,0.0,0.0,text);
  sprintf(text,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d",img.naxis1*sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1])/3600.0,img.naxis2*sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2])/3600.0,sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1]),sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2]));
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

// Computes apparent position
struct sat apparent_position(double mjd)
{
  struct sat s;
  double jd,rsun,rearth,rsat;
  double dx,dy,dz,dvx,dvy,dvz;
  xyz_t satpos,obspos,obsvel,satvel,sunpos;
  double ra,de;
  double mjd0=51544.5;

  // Sat ID
  s.Isat=Isat;

  // Get Julian Date
  jd=mjd+2400000.5;

  // Get positions
  obspos_xyz(mjd,&obspos,&obsvel);
  satpos_xyz(jd,&satpos,&satvel);
  sunpos_xyz(mjd,&sunpos);

  // Sat positions
  s.x=satpos.x;
  s.y=satpos.y;
  s.z=satpos.z;
  s.vx=satvel.x;
  s.vy=satvel.y;
  s.vz=satvel.z;

  // Sun position from satellite
  dx=-satpos.x+sunpos.x;  
  dy=-satpos.y+sunpos.y;
  dz=-satpos.z+sunpos.z;

  // Distances
  rsun=sqrt(dx*dx+dy*dy+dz*dz);
  rearth=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);

  // Angles
  s.psun=asin(696.0e3/rsun)*R2D;
  s.pearth=asin(6378.135/rearth)*R2D;
  s.p=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;

  // Visibility state
  if (s.p-s.pearth<-s.psun)
    strcpy(s.state,"eclipsed");
  else if (s.p-s.pearth>-s.psun && s.p-s.pearth<s.psun)
    strcpy(s.state,"umbra");
  else if (s.p-s.pearth>s.psun)
    strcpy(s.state,"sunlit");

  // Position differences
  dx=satpos.x-obspos.x;  
  dy=satpos.y-obspos.y;
  dz=satpos.z-obspos.z;
  dvx=satvel.x-obsvel.x;
  dvy=satvel.y-obsvel.y;
  dvz=satvel.z-obsvel.z;
  
  // Celestial position
  s.r=sqrt(dx*dx+dy*dy+dz*dz);
  s.v=(dvx*dx+dvy*dy+dvz*dz)/s.r;
  ra=modulo(atan2(dy,dx)*R2D,360.0);
  de=asin(dz/s.r)*R2D;

  // Precess
  precess(mjd,ra,de,mjd0,&s.ra,&s.de);

  // Phase
  s.phase=acos(((obspos.x-satpos.x)*(sunpos.x-satpos.x)+(obspos.y-satpos.y)*(sunpos.y-satpos.y)+(obspos.z-satpos.z)*(sunpos.z-satpos.z))/(rsun*s.r))*R2D;
	  
  // Magnitude
  if (strcmp(s.state,"sunlit")==0) 
    s.mag=STDMAG-15.0+5*log10(s.r)-2.5*log10(sin(s.phase*D2R)+(M_PI-s.phase*D2R)*cos(s.phase*D2R));
  else
    s.mag=15;

  /*     
  // Convert and project
  if (strcmp(m.orientation,"horizontal")==0) {
    equatorial2horizontal(mjd,s.ra,s.de,&s.azi,&s.alt);
    forward(s.azi,s.alt,&s.rx,&s.ry);
  } else if (strcmp(m.orientation,"equatorial")==0) {
    forward(s.ra,s.de,&s.rx,&s.ry);
  }
  */
  return s;
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

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
}

// Get a x and y from a RA and Decl
void forward(double ra0,double de0,double ra,double de,double *x,double *y)
{
  int i;
  char pcode[4]="TAN";
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;

  // Initialize Projection Parameters
  prj.flag=0;
  prj.r0=0.;
  for (i=0;i<10;prj.p[i++]=0.);

  // Initialize Reference Angles
  cel.ref[0]=ra0;
  cel.ref[1]=de0;
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;

  if (celset(pcode,&cel,&prj)) {
    printf("Error in Projection (celset)\n");
    return;
  } else {
    if (celfwd(pcode,ra,de,&cel,&phi,&theta,&prj,x,y)) {
      printf("Error in Projection (celfwd)\n");
      return;
    }
  }
  *x*=3600.;
  *y*=3600.;

  return;
}

// Get a RA and Decl from x and y
void reverse(double ra0,double de0,double x,double y,double *ra,double *de)
{
  int i;
  char pcode[4]="TAN";
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;

  x/=3600.;
  y/=3600.;

  // Initialize Projection Parameters
  prj.flag=0;
  prj.r0=0.;
  for (i=0;i<10;prj.p[i++]=0.);

  // Initialize Reference Angles
  cel.ref[0]=ra0;
  cel.ref[1]=de0;
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;

  if (celset(pcode,&cel,&prj)) {
    printf("Error in Projection (celset)\n");
    return;
  } else {
    if (celrev(pcode,x,y,&prj,&phi,&theta,&cel,ra,de)) {
      printf("Error in Projection (celrev)\n");
      return;
    }
  }
  return;
}
