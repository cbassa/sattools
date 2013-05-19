#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "cpgplot.h"
#include "cel.h"
#include "sgdp4h.h"

#define LIM 128
#define NMAX 256
#define MMAX 1024
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
  double alpha0,delta0,ra0,de0,azi0,alt0;
  double fov,mjd,gmst,w,wl,wb;
  float length;
  float minmag,maxmag,minrad,maxrad;
  char orientation[LIM],projection[4],observer[32];
  char nfd[LIM],starfile[LIM],tlefile[LIM],iodfile[LIM];
  char datadir[LIM],tledir[LIM];
  double lat,lng;
  double h,sra,sde,sazi,salt;
  float alt,timezone;
  float fw,fh;
  int level,grid,site_id;
  int leoflag,iodflag,iodpoint,visflag,planar,pssatno,psnr;
  float psrmin,psrmax,rvis;
} m;
struct sat {
  long Isat;
  char state[10];
  float mag,age;
  double jd;
  double dx,dy,dz;
  double x,y,z,vx,vy,vz;
  double rsun,rearth,h;
  double psun,pearth,p,phase;
  double r,v,ra,de;
  double azi,alt;
  double rx,ry;
};
struct star {
  double ra,de;
  float pmra,pmde;
  float mag;
};
struct observation {
  int ssn,site;
  char iod_line[LIM];
  double mjd,ra,de,azi,alt;
  double lng,lat;
  float elv;
  float dt,st,dr,sr,dx,dy,t;
  int flag;
};
int fgetline(FILE *,char *,int);
double modulo(double,double);
void reverse(double,double,double *,double *);
void forward(double,double,double *,double *);
void init_plot(char *,float,float);
void skymap_plot_renew(void);
double gmst(double);
double dgmst(double);
void skymap_plothorizontal_grid();
void skymap_plotequatorial_grid();
void skymap_plotconstellations(char *);
void equatorial2horizontal(double,double,double,double *,double *);
void horizontal2equatorial(double,double,double,double *,double *);
void skymap_plotstars(char *);
void obspos_xyz(double,xyz_t *,xyz_t *);
void sunpos_xyz(double,xyz_t *,double *,double *);
void skymap_plotsatellite(char *,int,double,double);
double date2mjd(int,int,double);
struct sat apparent_position(double);
long identify_satellite(char *,int,double,float,float);
int plot_skymap(void);
void rotate(int,float,float *,float *,float *);
int print_tle(char *,int); 
void mjd2date(double mjd,char *date);
void dec2sex(double x,char *s,int f,int len);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);
double nfd2mjd(char *date);
void nfd_now(char *s);
double sex2dec(char *s);
double doy2mjd(int year,double doy);
struct observation decode_iod_observation(char *iod_line);
void plot_iod(char *filename);
void get_site(int site_id);

void usage()
{
  
  printf("skymap t:c:i:R:D:hs:d:l:P:r:V:\n\n");
  printf("t    date/time (yyyy-mm-ddThh:mm:ss.sss) [default: now]\n");
  printf("c    TLE catalog file [default: classfd.tle]\n");
  printf("i    satellite ID (NORAD) [default: all]\n");
  printf("R    R.A.\n");
  printf("D    Decl.\n");
  printf("h    this help\n");
  printf("s    site (COSPAR)\n");
  printf("d    IOD observations\n");
  printf("l    trail length [default: 60s]\n");
  printf("P    planar search satellite ID\n");
  printf("r    planar search altitude\n");
  printf("V    altitude for visibility contours\n");

  return;
}

void init_skymap(void)
{
  int i;
  char *env;

  // Default Map parameters
  m.azi0=0;
  m.alt0=90.0;
  m.w=120.0;
  m.wl=180.0;
  m.wb=180.0;
  m.level=1;
  m.minmag=-2.0;
  m.maxmag=5.0;
  m.maxrad=2.0;
  m.minrad=0.02;
  strcpy(m.orientation,"horizontal");
  strcpy(m.starfile,"hip6mag.dat");
  strcpy(m.projection,"STG");

  m.lat=0.0;
  m.lng=0.0;
  m.alt=0.0;
  m.timezone=+0.0;
  m.grid=1;
  m.length=60.0;
  m.mjd=-1.0;
  m.leoflag=1;
  m.iodflag=0;
  m.visflag=0;
  m.planar=0;

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

void read_iod(char *filename,int iobs) 
{
  int i=0;
  char line[LIM];
  FILE *file;
  struct observation obs;

  file=fopen(filename,"r");
  // Read data
  while (fgets(line,LIM,file)!=NULL) {
    if (strlen(line)<10)
      continue;
    if (strstr(line,"#")==NULL) {
      obs=decode_iod_observation(line);
      if (i==iobs) {
	printf("%s\n",obs.iod_line);
	break;
      }
      i++;
    }
  }
  fclose(file);
 
  // Set parameters
  get_site(obs.site);
  m.mjd=obs.mjd;
  m.ra0=obs.ra;
  m.de0=obs.de;
  strcpy(m.orientation,"equatorial");
  m.level=6;

  return;
}


int main(int argc,char *argv[])
{
  int i,arg=0;

  // Redirect stderr
  freopen("/dev/null","w",stderr);

  init_skymap();

  // Decode options
  while ((arg=getopt(argc,argv,"t:c:i:R:D:hs:d:l:P:r:V:"))!=-1) {
    switch(arg) {

    case 't':
      strcpy(m.nfd,optarg);
      m.mjd=nfd2mjd(m.nfd);
      break;

    case 'c':
      strcpy(m.tlefile,optarg);
      break;

    case 'd':
      strcpy(m.iodfile,optarg);
      m.iodpoint=0;
      m.leoflag=0;
      read_iod(m.iodfile,m.iodpoint);
      m.iodflag=1;
      break;

    case 'l':
      m.length=atof(optarg);
      break;

    case 's':
      get_site(atoi(optarg));
      break;

    case 'i':
      Isatsel=atoi(optarg);
      m.leoflag=0;
      break;
      
    case 'P':
      m.planar=1;
      m.pssatno=atoi(optarg);
      m.psrmin=300;
      m.psrmax=1000;
      m.psnr=8;
      break;

    case 'r':
      m.psrmin=atof(optarg);
      m.psrmax=atof(optarg);
      m.psnr=1;
      break;

    case 'V':
      m.visflag=1;
      m.rvis=atof(optarg);
      break;

    case 'R':
      m.ra0=15.0*sex2dec(optarg);
      strcpy(m.orientation,"equatorial");
      m.level=5;
      break;

    case 'D':
      m.de0=sex2dec(optarg);
      strcpy(m.orientation,"equatorial");
      m.level=5;
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

  init_plot("/xs",10,0.75);
  
  plot_skymap();

  cpgend();

  fclose(stderr);

  return 0;
}

// Plot visibility contours
void plot_visibility(float h)
{
  int i,j,k,nx=300,ny=200,nc;
  float xmin,xmax,ymin,ymax;
  double rx,ry,azi,alt,ra,de;
  xyz_t obspos,obsvel,satpos,sunpos;
  double dx,dy,dz,r,ax,ay,az,d,dr,sra,sde;
  float rsun,rearth,psun,pearth,p,phase,mag;
  char state[10];
  float *cont,cmax;
  float tr[6];
  float c[]={0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};

  // Allocate
  cont=(float *) malloc(sizeof(float)*nx*ny);

  // Limits
  xmin=-1.5*m.w;
  xmax=1.5*m.w;
  ymin=-m.w;
  ymax=m.w;

  // Transformation matrix
  tr[2]=0.0;
  tr[1]=(xmax-xmin)/(float) nx;
  tr[0]=xmin-0.5*tr[1];
  tr[4]=0.0;
  tr[5]=(ymax-ymin)/(float) ny;
  tr[3]=ymin-0.5*tr[5];

  // Get observer and solar position
  obspos_xyz(m.mjd,&obspos,&obsvel);
  sunpos_xyz(m.mjd,&sunpos,&sra,&sde);

  for (i=0;i<nx;i++) {
    rx=xmin+(xmax-xmin)*(double) i/(double) (nx-1);
    for (j=0;j<ny;j++) {
      ry=ymin+(ymax-ymin)*(double) j/(double) (ny-1);
      reverse(rx,ry,&azi,&alt);
      
      // Skip low elevations
      if (alt<=0.0) {
	mag=15.0;
      } else {
	horizontal2equatorial(m.mjd,azi,alt,&ra,&de);
	
	// Compute unit vector
	ax=cos(ra*D2R)*cos(de*D2R);
	ay=sin(ra*D2R)*cos(de*D2R);
	az=sin(de*D2R);
	
	// Find distance
	for (k=0,d=h;k<20;k++) { 
	  dx=d*ax;
	  dy=d*ay;
	  dz=d*az;
	  satpos.x=obspos.x+dx;
	  satpos.y=obspos.y+dy;
	  satpos.z=obspos.z+dz;
	  r=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
	  dr=h+XKMPER-r;
	  if (dr<1.0)
	    break;
	  d+=dr;
	}
	
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
	p=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;
	p-=pearth;

	// Visibility
	if (p<-psun) {
	  strcpy(state,"eclipsed");
	  cpgsci(14);
	} else if (p>-psun && p<psun) {
	  strcpy(state,"umbra");
	  cpgsci(15);
	} else if (p>psun) {
	  strcpy(state,"sunlit");
	  cpgsci(7);
	}
	
	// Phase
	phase=acos(((obspos.x-satpos.x)*(sunpos.x-satpos.x)+(obspos.y-satpos.y)*(sunpos.y-satpos.y)+(obspos.z-satpos.z)*(sunpos.z-satpos.z))/(rsun*d))*R2D;
	
	// Magnitude
	if (strcmp(state,"sunlit")==0) 
	  mag=STDMAG-15.0+5*log10(d)-2.5*log10(sin(phase*D2R)+(M_PI-phase*D2R)*cos(phase*D2R));
	else
	  mag=15;
      }
      k=i+nx*j;
      cont[k]=mag;
    }
  }
  
  // Find maximum contour value
  for (i=0,j=0;i<nx*ny;i++) {
    if (cont[i]<15.0) {
      if (j==0 || cont[i]>cmax) cmax=cont[i];
      j++;
    }
  }
  nc=(int) cmax+1.0;

  // Plot contours
  cpgsci(7);
  cpgcont(cont,nx,ny,1,nx,1,ny,c,nc,tr);

  // Label contours
  cpgsch(0.8);
  for (i=0;i<nc;i++) {
    sprintf(state,"%.0f",c[i]);
    cpgconl(cont,nx,ny,1,nx,1,ny,c[i],tr,state,300,18);
  }
  cpgsch(1.0);
  cpgsci(1);

  return;
}

// Initialize plot
void init_plot(char *psfile,float width,float aspect)
{
  int i;

  // Initialize plot
  cpgopen(psfile);
  cpgslw(2);
  cpgpap(width,aspect);

  skymap_plot_renew();

  return;
}

// Add to schedule
void schedule(char *nfd,double ra,double de)
{
  FILE *file;
  char sra[16],sde[16];

  // Compute strings
  dec2sex(ra/15.0,sra,0,5);
  dec2sex(de,sde,0,4);

  printf("%s %s %s\n",nfd,sra,sde);

  // Open file
  file=fopen("schedule.txt","a");
  if (file==NULL) {
    printf("Failed to create schedule.txt\n");
    return;
  }
  fprintf(file,"%s %s %s\n",nfd,sra,sde);
  fclose(file);

  return;
}

// Initialize plot
void skymap_plot_renew(void)
{
  char filename[LIM];

  // Size limit
  if (m.w>120.0)
    m.w=120.0;

  if (m.level==1) {
    m.w=120;
    m.minmag=-2.0;
    m.maxmag=5.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/hip6mag.dat",m.datadir);
  }
  if (m.level==2) {
    m.w=90;
    m.minmag=-1.5;
    m.maxmag=5.5;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/hip6mag.dat",m.datadir);
  }
  if (m.level==3) {
    m.w=60;
    m.minmag=-1.0;
    m.maxmag=6.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/hip6mag.dat",m.datadir);
  }
  if (m.level==4) {
    m.w=30;
    m.minmag=-0.5;
    m.maxmag=6.5;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/tyc8mag.dat",m.datadir);
  }
  if (m.level==5) {
    m.w=20;
    m.minmag=0.0;
    m.maxmag=7.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/tyc8mag.dat",m.datadir);
  }
  if (m.level==6) {
    m.w=10;
    m.minmag=1.0;
    m.maxmag=8.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/tyc8mag.dat",m.datadir);
  }
  if (m.level==7) {
    m.w=5;
    m.minmag=2.0;
    m.maxmag=9.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/tyc10mag.dat",m.datadir);
  }
  if (m.level==8) {
    m.w=2;
    m.minmag=3.0;
    m.maxmag=10.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/tyc10mag.dat",m.datadir);
  }
  if (m.level==9) {
    m.w=1;
    m.minmag=3.0;
    m.maxmag=10.0;
    m.maxrad=2.0;
    m.minrad=0.02;
    sprintf(m.starfile,"%s/data/tyc12mag.dat",m.datadir);
  }
  /*
  // Star files
  if (m.w>90.0 && m.w<=120.0) {
    strcpy(m.starfile,"data/hip6mag.dat");
    m.minmag=-2.0;
    m.maxmag=4.5;
  } else if (m.w>60.0 && m.w<=90.0) {
    strcpy(m.starfile,"data/hip6mag.dat");
    m.minmag=-2.0;
    m.maxmag=5.0;
  } else if (m.w>20.0 && m.w<=60.0) {
    strcpy(m.starfile,"data/hip6mag.dat");
    m.minmag=0.0;
    m.maxmag=6.0;
  } else if (m.w<=20.0) {
    strcpy(m.starfile,"data/tyc8mag.dat");
    m.minmag=2.0;
    m.maxmag=8.0;
  }
  */


  return;
}

// Get a x and y from an AZI, ALT
void forward(double alpha,double delta,double *x,double *y)
{
  int i;
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;

  // Initialize Projection Parameters
  prj.flag=0;
  prj.r0=0.;
  for (i=0;i<10;prj.p[i++]=0.);

  // Initialize Reference Angles
  if (strcmp(m.orientation,"horizontal")==0) {
    cel.ref[0]=m.azi0;
    cel.ref[1]=m.alt0;
  } else if (strcmp(m.orientation,"equatorial")==0) {
    cel.ref[0]=m.ra0;
    cel.ref[1]=m.de0;
  }
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;

  if (celset(m.projection,&cel,&prj)) {
    printf("Error in Projection (celset)\n");
    return;
  } else {
    if (celfwd(m.projection,alpha,delta,&cel,&phi,&theta,&prj,x,y)) {
      printf("Error in Projection (celfwd)\n");
      return;
    }
  }

  // Flip equatorial axis
  if (strcmp(m.orientation,"equatorial")==0)
    *x*=-1;

  return;
}

// Get an AZI, ALT from x and y
void reverse(double x,double y,double *alpha,double *delta)
{
  int i;
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;

  // Flip equatorial axis
  if (strcmp(m.orientation,"equatorial")==0)
    x*=-1;

  // Initialize Projection Parameters
  prj.flag=0;
  prj.r0=0.;
  for (i=0;i<10;prj.p[i++]=0.);

  // Initialize Reference Angless
  if (strcmp(m.orientation,"horizontal")==0) {
    cel.ref[0]=m.azi0;
    cel.ref[1]=m.alt0;
  } else if (strcmp(m.orientation,"equatorial")==0) {
    cel.ref[0]=m.ra0;
    cel.ref[1]=m.de0;
  }
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;

  if (celset(m.projection,&cel,&prj)) {
    printf("Error in Projection (celset)\n");
    return;
  } else {
    if (celrev(m.projection,x,y,&prj,&phi,&theta,&cel,alpha,delta)) {
      printf("Error in Projection (celrev)\n");
      return;
    }
  }
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

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
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

// Plot field of view
void skymap_plot_fov()
{
  int i,j,n=50;
  double rx,ry,azi,alt;
  float x,y;
  //  double azi0[]={44.0,7.0,30.0,44.0,44.0};
  //  double alt0[]={16.0,16.0,60.0,60.0,16.0};
  double azi0[]={45.0,7.0,7.0,45.0,45.0};
  double alt0[]={20.0,20.0,50.0,50.0,20.0};
  double azi1[]={150.0,140.0,140.0,150.0,150.0};
  double alt1[]={26.0,26.0,31.0,31.0,26.0};

  /*
  for (i=0;i<sizeof(alt0)/sizeof(alt0[0])-1;i++) {
    for (j=0;j<n;j++) {
      azi=azi0[i]+(azi0[i+1]-azi0[i])*(float) j/(float) (n-1);
      alt=alt0[i]+(alt0[i+1]-alt0[i])*(float) j/(float) (n-1);
      azi+=180.0;
      forward(azi,alt,&rx,&ry);
      x=(float) rx;
      y=(float) ry;
      if (j==0)
	cpgmove(x,y);
      else
	cpgdraw(x,y);
    }
  }
  */
  for (i=0;i<sizeof(alt1)/sizeof(alt1[0])-1;i++) {
    for (j=0;j<n;j++) {
      azi=azi1[i]+(azi1[i+1]-azi1[i])*(float) j/(float) (n-1);
      alt=alt1[i]+(alt1[i+1]-alt1[i])*(float) j/(float) (n-1);
      azi+=180.0;
      forward(azi,alt,&rx,&ry);
      x=(float) rx;
      y=(float) ry;
      if (j==0)
	cpgmove(x,y);
      else
	cpgdraw(x,y);
    }
  }

  
  return;
}

// Plot an horizontal grid
void skymap_plothorizontal_grid()
{
  int i,j;
  double rx,ry,azi,alt;
  float x,y;
  float sch;
  int sci,sls,status;
  FILE *file;
  char line[LIM],filename[LIM];

  // Get setup
  cpgqch(&sch);
  cpgqls(&sls);
  cpgqci(&sci);

  // Plot grid
  cpgsci(15);
  cpgsls(2);

  // Altitudes
  if (m.grid==1) {
    for (alt=0.0;alt<=80.0;alt+=20.0) {
      for (i=0;i<NMAX;i++) {
	azi=360.0*(double) i/(double) (NMAX-1);
	
	forward(azi,alt,&rx,&ry);
	x=(float) rx;
	y=(float) ry;
	if (i==0) cpgmove(x,y);
	if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
	  cpgdraw(x,y);
	else
	  cpgmove(x,y);
      }
    }
    
    // Azimuths
    for (azi=0.0;azi<360.0;azi+=30.0) {
      for (i=0;i<NMAX;i++) {
	alt=0.0+80.0*(double) i/(double) (NMAX-1);
	
	forward(azi,alt,&rx,&ry);
	x=(float) rx;
	y=(float) ry;
        if (i==0) cpgmove(x,y);
	if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
	  cpgdraw(x,y);
	else
	  cpgmove(x,y);
      }
    }
  }
  cpgsci(1);
  cpgsls(1);

  // Plot horizon
  for (i=0;i<NMAX;i++) {
    azi=360.0*(float) i/(float) (NMAX-1);
    alt=0.0;
    
    forward(azi,alt,&rx,&ry);
    x=(float) rx;
    y=(float) ry;
    if (i==0) cpgmove(x,y);
    if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
      cpgdraw(x,y);
    else
      cpgmove(x,y);
  }

  // Use real horizon
  sprintf(filename,"%s/data/%04dhorizon.txt",m.datadir,m.site_id);
  file=fopen(filename,"r");
  if (file!=NULL) {
    i=0;
    while (fgetline(file,line,LIM)>0) {
      if (strlen(line)<2) {
	i=0;
	continue;
      }
      status=sscanf(line,"%lf %lf",&azi,&alt);
      forward(azi+180.0,alt,&rx,&ry);
      x=(float) rx;
      y=(float) ry;
      if (i==0) cpgmove(x,y);
      if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
	cpgdraw(x,y);
      else
	cpgmove(x,y);
      i++;
    }
    fclose(file);
  }

  return;
}

// Plot an equatorial grid
void skymap_plotequatorial_grid()
{
  int i,j;
  double rx,ry,ra,de;
  float x,y;
  float sch;
  int sci,sls;

  // Get setup
  cpgqch(&sch);
  cpgqls(&sls);
  cpgqci(&sci);

  // Plot grid
  cpgsci(15);
  cpgsls(2);

  // Declinations
  if (m.grid==1) {
    for (de=-80.0;de<=80.0;de+=20.0) {
      for (i=0;i<NMAX;i++) {
	ra=360.0*(double) i/(double) (NMAX-1);
	
	forward(ra,de,&rx,&ry);
	x=(float) rx;
	y=(float) ry;
	if (i==0) cpgmove(x,y);
	if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
	  cpgdraw(x,y);
	else
	  cpgmove(x,y);
      }
    }
    
    // Right ascensions
    for (ra=0.0;ra<360.0;ra+=30.0) {
      for (i=0;i<NMAX;i++) {
	de=-80.0+160.0*(double) i/(double) (NMAX-1);
	
	forward(ra,de,&rx,&ry);
	x=(float) rx;
	y=(float) ry;
        if (i==0) cpgmove(x,y);
	if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
	  cpgdraw(x,y);
	else
	  cpgmove(x,y);
      }
    }
  }
  cpgsci(1);
  cpgsls(1);

  // Plot equatorial
  for (i=0;i<NMAX;i++) {
    ra=360.0*(float) i/(float) (NMAX-1);
    de=0.0;
    
    forward(ra,de,&rx,&ry);
    x=(float) rx;
    y=(float) ry;
    if (i==0) cpgmove(x,y);
    if (fabs(x)<=1.5*m.w && fabs(y)<=m.w)
      cpgdraw(x,y);
    else
      cpgmove(x,y);
  }
  

  return;
}

// Plot constellation lines
void skymap_plotconstellations(char *filename)
{
  int i,flag;
  unsigned long tyc;
  double rx,ry,azi,alt,alt1;
  double ra,de,ra0,de0,mjd0=51544.5;
  //  double mjd0=48348.3125; // J1991.25
  float x,y,x1,y1;
  char line[LIM];
  FILE *file;

  cpgsci(4);
  // Loop over file
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("Const file not found\n");
    exit(1);
  }
    
  for (i=0;fgetline(file,line,LIM)>0;i++) {
    if (strchr(line,'#')!=NULL) continue;
    sscanf(line,"%lu %i %lf %lf\n",&tyc,&flag,&ra0,&de0);
    precess(mjd0,ra0,de0,m.mjd,&ra,&de);
    if (strcmp(m.orientation,"horizontal")==0) {
      equatorial2horizontal(m.mjd,ra,de,&azi,&alt);
      forward(azi,alt,&rx,&ry);
    } else if (strcmp(m.orientation,"equatorial")==0) {
      forward(ra,de,&rx,&ry);
    }
    x=(float) rx;
    y=(float) ry;
    if (i==0) cpgmove(x,y);
    //    if (fabs(x)<=m.w && fabs(y)<=m.w) {
      if (flag==0) cpgmove(x,y);
      if (flag==1) cpgdraw(x,y);
      //    } else
      //      cpgmove(x,y);
  }
  fclose(file);
  cpgsci(1);

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


// Plot the stars
void skymap_plotstars(char *filename)
{
  int i,j;
  unsigned long hip;
  double rx,ry,azi,alt;
  double ra,de,ra0,de0,mjd0=51544.5;
  //  double mjd0=48348.3125; // J1991.25
  float x,y,vmag;
  float rad;
  char line[LIM];
  FILE *file;
  struct star s;

  /*
  // 20130118 Testing of the binary catalog also used in plotfits
  file=fopen("software/sattools/data/tycho2.dat","rb");
  while (!feof(file)) {
    fread(&s,sizeof(struct star),1,file);
    vmag=s.mag;
    ra0=s.ra;
    de0=s.de;
    if (vmag<=m.maxmag) {
      precess(mjd0,ra0,de0,m.mjd,&ra,&de);
      if (strcmp(m.orientation,"horizontal")==0) {
	equatorial2horizontal(m.mjd,ra,de,&azi,&alt);
	forward(azi,alt,&rx,&ry);
      } else if (strcmp(m.orientation,"equatorial")==0) {
	forward(ra,de,&rx,&ry);
      }
      x=(float) rx;
      y=(float) ry;
      
      // Star size
      rad=m.maxrad+(m.minrad-m.maxrad)*(vmag-m.minmag)/(m.maxmag-m.minmag);
      rad*=m.w/90.0;
      cpgsci(0);
      cpgcirc(x,y,1.3*rad);
      cpgsci(1);
      cpgcirc(x,y,rad);
    }
  }
  fclose(file);
  */
    
  // Loop over file
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("Star file not found\n");
    exit(1);
  }
  while (fgetline(file,line,LIM)>0) {
    // Failed for Tycho files
    //    sscanf(line,"%i %lf %lf %f\n",&hip,&ra,&de,&vmag);
    // Skipping star ID
    sscanf(line+13,"%lf %lf %f\n",&ra0,&de0,&vmag);

    if (vmag<=m.maxmag) {
      precess(mjd0,ra0,de0,m.mjd,&ra,&de);
      if (strcmp(m.orientation,"horizontal")==0) {
	equatorial2horizontal(m.mjd,ra,de,&azi,&alt);
	forward(azi,alt,&rx,&ry);
      } else if (strcmp(m.orientation,"equatorial")==0) {
	forward(ra,de,&rx,&ry);
      }
      x=(float) rx;
      y=(float) ry;
      
      if (fabs(rx)<0.02 && fabs(ry)<0.02)
	printf("%lf %lf %lf %lf %f\n",ra0,de0,ra,de,vmag);

      // Star size
      rad=m.maxrad+(m.minrad-m.maxrad)*(vmag-m.minmag)/(m.maxmag-m.minmag);
      rad*=m.w/90.0;
      cpgsci(0);
      cpgcirc(x,y,1.3*rad);
      cpgsci(1);
      cpgcirc(x,y,rad);
    }
  }
  fclose(file);
  
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

// Plot satellite track
void skymap_plotsatellite(char *filename,int satno,double mjd0,double dt)
{
  orbit_t orb;
  struct sat s;
  int imode,flag,fflag,i;
  FILE *fp=NULL;
  xyz_t satpos,obspos,satvel,sunpos;
  double mjd,jd,dx,dy,dz;
  double rx,ry,ra,de,azi,alt,r,t;
  float x,y;
  char norad[7],satname[30],date[24];;
  float isch;
  float rsun,rearth,psun,pearth,p;
  int priority[]={24680,28888,15071,26934,37348,5678,5679,5680,5681,5682,8818,8835,8836,8884,10502,10529,10544,10594,11720,11731,11732,11745,13791,13844,13845,13874};

  // Open TLE file
  fp=fopen(filename,"rb");
  if (fp==NULL)
    fatal_error("File open failed for reading %s\n",filename);

  // Read TLEs
  while (read_twoline(fp,satno,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);

    if (m.leoflag==1 && orb.rev<3)
      continue;
    if (m.leoflag==2 && orb.rev>=3)
      continue;

    sprintf(norad," %ld",Isat);

    if (imode==SGDP4_ERROR)
      continue;

    for (flag=0,fflag=0,t=0.0;t<dt;t+=1.0) {
      mjd=mjd0+t/86400.0;

      // Compute apparent position
      s=apparent_position(mjd);

      // Convert to float
      x=(float) s.rx;
      y=(float) s.ry;

      // Visibility
      if (s.p-s.pearth<-s.psun) {
	cpgsci(14);
      } else if (s.p-s.pearth>-s.psun && s.p-s.pearth<s.psun) {
	cpgsci(15);
      } else if (s.p-s.pearth>s.psun) {
	for (i=0;i<sizeof(priority)/sizeof(priority[0]);i++) {
	  if (Isat==priority[i]) {
	    cpgsci(2);
	    break;
	  } else
	    cpgsci(7);
	}
      }

      // In field of view
      if (fabs(x)<m.fw && fabs(y)<m.fh && fflag==0) {
	mjd2date(mjd,date);
	printf("%.19s %s %6.1f\n",date,norad,s.mag);
	fflag=1;
      }

      // Plot satellites
      if (flag==0) {
	if (s.age<25) 
	  cpgpt1(x,y,17);
	else if (s.age<50)
	  cpgpt1(x,y,4);
	else 
	  cpgpt1(x,y,6);
	cpgsch(0.6);

	// Print name if in viewport
	if (fabs(x)<1.5*m.w && fabs(y)<m.w && x<1.32*m.w && y<0.96*m.w)
	  cpgtext(x,y,norad);
	cpgsch(isch);
	cpgmove(x,y);
	flag=1;
      } else {
	cpgdraw(x,y);
      }
    }
  }
  fclose(fp);
  cpgsci(1); 

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


// Computes apparent position
struct sat apparent_position(double mjd)
{
  struct sat s;
  double jd,rsun,rearth,rsat;
  double dx,dy,dz,dvx,dvy,dvz;
  xyz_t satpos,obspos,obsvel,satvel,sunpos;
  double sra,sde;

  // Sat ID
  s.Isat=Isat;

  // Get Julian Date
  jd=mjd+2400000.5;

  // Get positions
  obspos_xyz(mjd,&obspos,&obsvel);
  satpos_xyz(jd,&satpos,&satvel);
  sunpos_xyz(mjd,&sunpos,&sra,&sde);

  // Age
  s.age=jd-SGDP4_jd0;

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
  s.h=rearth-XKMPER;
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
  s.ra=modulo(atan2(dy,dx)*R2D,360.0);
  s.de=asin(dz/s.r)*R2D;

  // Phase
  s.phase=acos(((obspos.x-satpos.x)*(sunpos.x-satpos.x)+(obspos.y-satpos.y)*(sunpos.y-satpos.y)+(obspos.z-satpos.z)*(sunpos.z-satpos.z))/(rsun*s.r))*R2D;
	  
  // Magnitude
  if (strcmp(s.state,"sunlit")==0) 
    s.mag=STDMAG-15.0+5*log10(s.r)-2.5*log10(sin(s.phase*D2R)+(M_PI-s.phase*D2R)*cos(s.phase*D2R));
  else
    s.mag=15;

     
  // Convert and project
  if (strcmp(m.orientation,"horizontal")==0) {
    equatorial2horizontal(mjd,s.ra,s.de,&s.azi,&s.alt);
    forward(s.azi,s.alt,&s.rx,&s.ry);
  } else if (strcmp(m.orientation,"equatorial")==0) {
    forward(s.ra,s.de,&s.rx,&s.ry);
  }

  return s;
}

// Planar search
void planar_search(char *filename,int satno,float rmin,float rmax,int nr)
{
  int i,j,imode;
  FILE *fp;
  orbit_t orb;
  kep_t K;
  int withvel,rv;
  double tsince,radius,jd;
  double st,ct,sn,cn,si,ci,t;
  xyz_t satpos,obspos,obsvel,sunpos;
  double r,ra,de,dx,dy,dz,rsun,rearth,psun,pearth,p,azi,alt,rx,ry,rx0,ry0,ra0,de0;
  double sra,sde;
  float phase,mag,mmin;
  char state[10];

  // Open TLE file
  fp=fopen(filename,"rb");
  if (fp==NULL)
    fatal_error("File open failed for reading %s\n",filename);

  // Read TLEs
  while (read_twoline(fp,satno,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);

    if (imode==SGDP4_ERROR)
      continue;

  }
  fclose(fp);

  // Get Julian Date
  jd=m.mjd+2400000.5;

  // Get kepler
  tsince=1440.0*(jd-SGDP4_jd0);
  rv=sgdp4(tsince,1,&K);

  // Angles
  sn=sin(K.ascn);
  cn=cos(K.ascn);
  si=sin(K.eqinc);
  ci=cos(K.eqinc);

  // Loop over radii
  for (j=0;j<nr;j++) {
    if (nr>1)
      radius=rmin+(rmax-rmin)*(float) j/(float) (nr-1);
    else
      radius=rmin;

    // Loop over angles
    for (i=0,mmin=15.0;i<MMAX;i++) {
      t=2.0*M_PI*(double) i/(double) (MMAX-1);
      st=sin(K.theta+t);
      ct=cos(K.theta+t);
      
      satpos.x=-sn*ci*st+cn*ct;
      satpos.y=cn*ci*st+sn*ct;
      satpos.z=si*st;
      satpos.x*=(radius+XKMPER);
      satpos.y*=(radius+XKMPER);
      satpos.z*=(radius+XKMPER);
      
      // Get positions
      obspos_xyz(m.mjd,&obspos,&obsvel);
      sunpos_xyz(m.mjd,&sunpos,&sra,&sde);
      
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
      p=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;
      p-=pearth;
      
      // Position differences
      dx=satpos.x-obspos.x;  
      dy=satpos.y-obspos.y;
      dz=satpos.z-obspos.z;
      
      // Celestial position
      r=sqrt(dx*dx+dy*dy+dz*dz);
      ra=atan2(dy,dx)*R2D;
      de=asin(dz/r)*R2D;
      
      // Convert and project
      if (strcmp(m.orientation,"horizontal")==0) {
	equatorial2horizontal(m.mjd,ra,de,&azi,&alt);
	forward(azi,alt,&rx,&ry);
      } else if (strcmp(m.orientation,"equatorial")==0) {
	forward(ra,de,&rx,&ry);
      }
      
      // Visibility
      if (p<-psun) {
	strcpy(state,"eclipsed");
	cpgsci(14);
      } else if (p>-psun && p<psun) {
	strcpy(state,"umbra");
	cpgsci(15);
      } else if (p>psun) {
	strcpy(state,"sunlit");
	cpgsci(7);
      } 

      // Phase
      phase=acos(((obspos.x-satpos.x)*(sunpos.x-satpos.x)+(obspos.y-satpos.y)*(sunpos.y-satpos.y)+(obspos.z-satpos.z)*(sunpos.z-satpos.z))/(rsun*r))*R2D;
      
      // Magnitude
      if (strcmp(state,"sunlit")==0) 
	mag=STDMAG-15.0+5*log10(r)-2.5*log10(sin(phase*D2R)+(M_PI-phase*D2R)*cos(phase*D2R));
      else
	mag=15;
      
      if (mag<mmin) {
	rx0=rx;
	ry0=ry;
	mmin=mag;
      }

      if (i==0)
	cpgmove((float) rx,(float) ry);
      else
	cpgdraw((float) rx,(float) ry);
      cpgsci(1);
    }
    
    // Plot brightest point
    cpgsch(1.0);
    cpgsci(7);
    cpgpt1((float) rx0,(float) ry0,4);
    cpgsci(1);
  }

  return;
}

// Print TLE
int print_tle(char *filename,int satno) 
{
  FILE *fp;
  struct sat s;
  orbit_t orb;
  char line[LIM],pline[LIM];
  int flag=0;

  // Open TLE file
  fp=fopen(filename,"rb");
  if (fp==NULL)
    fatal_error("File open failed for reading %s\n",filename);

  // Read TLEs
  while (fgetline(fp,line,LIM)>0) {
    sscanf(line+2,"%ld",&Isat);
    if (line[0]=='1' && (long) Isat==satno) {
      printf("\n%s%s",pline,line);
      fgetline(fp,line,LIM);
      printf("%s\n",line);
      flag=1;
    }

    strcpy(pline,line);
  }
  fclose(fp);

  return flag;
}

// Identify satellite
long identify_satellite(char *filename,int satno,double mjd,float rx,float ry)
{
  long Isatmin=0;
  int imode;
  FILE *fp;
  struct sat s,smin;
  orbit_t orb;
  float dr,drmin,agemin;
  char line[LIM],pline[LIM];
  char sra[16],sde[16];

  // Open TLE file
  fp=fopen(filename,"rb");
  if (fp==NULL)
    fatal_error("File open failed for reading %s\n",filename);

  // Read TLEs
  while (read_twoline(fp,satno,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);

    if (imode==SGDP4_ERROR)
      continue;

    // Compute apparent position
    s=apparent_position(mjd);

    // Offset
    dr=sqrt(pow(rx-s.rx,2)+pow(ry-s.ry,2));

    if (Isatmin==0 || dr<drmin) {
      Isatmin=Isat;
      drmin=dr;
      smin=s;
    }
  }
  fclose(fp);

  // Print TLE
  print_tle(filename,Isatmin);
  printf("Age: %.1f d\n\n",smin.age);
  printf("x: %+10.2lf km; vx: %+6.3f km/s\ny: %+10.2lf km; vy: %+6.3f km/s\nz: %+10.2lf km; vz: %+6.3f km/s\nr: %10.2lf km; v:  %6.3f km/s\nh: %10.2lf km\n\n",smin.x,smin.vx,smin.y,smin.vy,smin.z,smin.vz,smin.r,smin.v,smin.h);
  dec2sex(smin.ra/15.0,sra,0,5);
  dec2sex(smin.de,sde,0,4);

  printf("R.A.: %s  Decl.: %s\n",sra,sde);
  printf("Azi.: %.1f Alt.: %.1f\n\n",modulo(smin.azi-180.0,360.0),smin.alt);

  printf("Phase: %.2f\nMagnitude: %.2f\n",smin.phase,smin.mag);
  
  return Isatmin;
}

void skymap_plotsun(void)
{
  double rx,ry;

  if (strcmp(m.orientation,"horizontal")==0) 
    forward(m.sazi,m.salt,&rx,&ry);
  else if (strcmp(m.orientation,"equatorial")==0) 
    forward(m.sra,m.sde,&rx,&ry);

  cpgsci(7);
  if (m.w>60.0)
    cpgcirc((float) rx,(float) ry,2.0);
  else
    cpgcirc((float) rx,(float) ry,0.25);
  cpgsci(1);

  return;
}

// plot skymap
int plot_skymap(void)
{
  int redraw=1,fov=2,status;
  float x,y;
  char c,text[256],sra[16],sde[16],filename[LIM];
  double ra,de,azi,alt,rx,ry;
  xyz_t sunpos;
  float focallength[]={28,35,50,100,200,300};

  for (;;) {
    if (redraw>0) {
      // Get present mjd
      if (m.mjd<0.0) {
	nfd_now(m.nfd);
	m.mjd=nfd2mjd(m.nfd);
      }
      
      // Get locations
      if (strcmp(m.orientation,"horizontal")==0) 
	horizontal2equatorial(m.mjd,m.azi0,m.alt0,&m.ra0,&m.de0);
      else if (strcmp(m.orientation,"equatorial")==0) 
	equatorial2horizontal(m.mjd,m.ra0,m.de0,&m.azi0,&m.alt0);
      
      // Get sun position
      sunpos_xyz(m.mjd,&sunpos,&m.sra,&m.sde);
      equatorial2horizontal(m.mjd,m.sra,m.sde,&m.sazi,&m.salt);

      cpgscr(0,0.0,0.0,0.0);
      cpgeras();

      // Create window
      cpgsvp(0.01,0.99,0.01,0.99);
      cpgwnad(-1.5*m.w,1.5*m.w,-m.w,m.w);

      // Set background
      if (m.salt>0.0) 
	cpgscr(0,0.0,0.0,0.4);
      else if (m.salt>-6.0)
	cpgscr(0,0.0,0.0,0.3);
      else if (m.salt>-12.0)
	cpgscr(0,0.0,0.0,0.2);
      else if (m.salt>-18.0)
	cpgscr(0,0.0,0.0,0.1);
      else
	cpgscr(0,0.0,0.0,0.0);
      cpgsci(0);
      cpgrect(-1.5*m.w,1.5*m.w,-m.w,m.w);
      cpgsci(1);

      cpgbox("BC",0.,0,"BC",0.,0);
      cpgpt1(0.0,0.0,2);

      // Plot field-of-view
      if (fov>=0) {
	cpgsfs(2);

	m.fw=atan(0.5*6.3265/focallength[fov])*R2D;
	m.fh=atan(0.5*4.6389/focallength[fov])*R2D;
	//m.fw=atan(0.5*22.3/focallength[fov])*R2D;
	//m.fh=atan(0.5*14.9/focallength[fov])*R2D;

	cpgrect(-m.fw,m.fw,-m.fh,m.fh);
	cpgsfs(1);
      }

      
      // Top left string
      cpgsch(0.8);
      mjd2date(m.mjd,m.nfd);
      sprintf(text,"%s UTC; %s (%04d) [%+.4f\\u\\(2218)\\d, %+.4f\\u\\(2218)\\d, %.0fm]",m.nfd,m.observer,m.site_id,m.lat,m.lng,m.alt*1000.0);
      cpgmtxt("T",0.6,0.0,0.0,text);
      
      // Top right string
      if (m.planar==0) {
	if (Isatsel==0) {
	  if (m.leoflag==-1)
	    sprintf(text,"None");
	  else if (m.leoflag==0)
	    sprintf(text,"All");
	  else if (m.leoflag==1)
	    sprintf(text,"LEO");
	  else if (m.leoflag==2)
	    sprintf(text,"HEO/GEO");
	} else if (Isatsel>0) {
	  sprintf(text,"%05d",(int) Isatsel);
	} else {
	  strcpy(text,"");
	}
      } else {
	if (Isatsel==0) {
	  if (m.leoflag==0)
	    sprintf(text,"Planar search: %05d; All",m.pssatno);
	  else if (m.leoflag==1)
	    sprintf(text,"Planar search: %05d; LEO",m.pssatno);
	  else if (m.leoflag==2)
	    sprintf(text,"Planar search: %05d; HEO/GEO",m.pssatno);
	} else if (Isatsel>0) {
	  sprintf(text,"Planar search: %05d; %05d",m.pssatno,(int) Isatsel);
	} else {
	  sprintf(text,"Planar search: %05d",m.pssatno);
	}
      }
      cpgmtxt("T",0.6,1.0,1.0,text);

      // Bottom string
      dec2sex(m.ra0/15.0,sra,0,5);
      dec2sex(m.de0,sde,0,4);
      sprintf(text,"R: %s; D: %s; A: %.1f; E: %.1f; S: %.1fx%.1f deg; L: %d; O: %s; m < %.1f; f: %.0f mm; l: %.0f s",sra,sde,modulo(m.azi0-180.0,360.0),m.alt0,3.0*m.w,2.0*m.w,m.level,m.orientation,m.maxmag,focallength[fov],m.length);
      cpgmtxt("B",1.0,0.0,0.0,text);
      cpgsch(1.0);

      // Plot everything
      if (strcmp(m.orientation,"horizontal")==0) {
       	skymap_plothorizontal_grid();
	horizontal2equatorial(m.mjd,m.azi0,m.alt0,&m.ra0,&m.de0);
	//	skymap_plot_fov();
      } else if (strcmp(m.orientation,"equatorial")==0) {
	skymap_plotequatorial_grid();
	equatorial2horizontal(m.mjd,m.ra0,m.de0,&m.azi0,&m.alt0);
      }
      sprintf(filename,"%s/data/constfig.dat",m.datadir);
      skymap_plotconstellations(filename);
      skymap_plotstars(m.starfile);
      
      if (Isatsel>=0 && m.leoflag>=0)
	skymap_plotsatellite(m.tlefile,Isatsel,m.mjd,m.length);

      skymap_plotsun();
    }
    // Reset redraw
    redraw=0;

    // Plot planar search
    if (m.planar==1) {
      if (m.pssatno==0) 
	printf("Please select a satellite.\n");
      else
	planar_search(m.tlefile,m.pssatno,m.psrmin,m.psrmax,m.psnr);
    }
    
    // Plot IOD points
    if (m.iodflag==1)
      plot_iod(m.iodfile);

    // Plot visibility
    if (m.visflag==1 && strcmp(m.orientation,"horizontal")==0)
      plot_visibility(m.rvis);

    // Get time
    cpgcurs(&x,&y,&c);
    
    // Help
    if (c=='h' || c=='H') {
      printf("q   quit\n");
      printf("i   Identify satellite\n");
      printf("r   Reset satellite selection/real time\n");
      printf("f   Select satellite\n");
      printf("l   Set integration length\n");
      printf("m   Measure cursor RA/Dec, Alt/Azi\n");
      printf("g   Toggle grid (on/off)\n");
      printf("o   Toggle orientation (horizontal/equatorial)\n");
      printf("c   Center on cursor\n");
      printf("z   Center on zenith\n");
      printf("n   Center on North\n");
      printf("s   Center on South\n");
      printf("e   Center on East\n");
      printf("w   Center on West\n");
      printf("1-9 Zoom level\n");
      printf("+   Zoom in one level\n");
      printf("-   Zoom out one level\n");
      printf(".   Increase time by 1 step\n");
      printf(",   Decrease time by 1 step\n");
      printf(">   Increase step size\n");
      printf("<   Decrease step size\n");
      printf("P   Toggle planar search\n");
      printf("R   Read catalog\n");
      printf("L   Toggle satellite selection (All, LEO, HEO/GEO, none)\n");
      printf("v   Toggle visibility contours\n");
      printf("F   Toggle focal length\n");
      printf("TAB Cycle IOD observations\n");
      printf("S   Save position/time to schedule\n");
    }

    // Cycle IOD points
    if (c=='\t') {
      m.iodpoint++;
      read_iod(m.iodfile,m.iodpoint);
      redraw=1;
    }

    // Toggle planar search
    if (c=='P') {
      if (m.planar==0) {
	printf("Provide altitude range (km): [min max num] ");
	status=scanf("%f %f %d",&m.psrmin,&m.psrmax,&m.psnr);
	m.pssatno=Isatsel;
	m.planar=1;
      } else if (m.planar==1) {
	m.pssatno=0;
	m.planar=0;
      }
      redraw=1;
    }

      

    // Toggle visibility contours
    if (c=='v') {
      if (m.visflag==0) {
	printf("Provide altitude (km): ");
	status=scanf("%f",&m.rvis);
	m.visflag=1;
      } else if (m.visflag==1) {
	m.visflag=0;
      }
      redraw=1;
    }

    // Identify
    if (c=='i') 
      identify_satellite(m.tlefile,Isatsel,m.mjd,x,y);
    

    // Read catalog
    if (c=='R') {
      printf("TLE catalog name: ");
      status=scanf("%s",m.tlefile);
      redraw=1;
    }

    // Increase/decrease time
    if (c=='.') {
      m.mjd+=m.length/86400.0;
      redraw=1;
    }      
    if (c==',') {
      m.mjd-=m.length/86400.0;
      redraw=1;
    }
    
    // Increase/decrease step
    if (c=='>') {
      m.length*=2.0;
      redraw=2;
    }
    if (c=='<') {
      m.length/=2.0;
      redraw=2;
    }
    
    // Reset
    if (c=='r') {
      Isatsel=0;
      m.length=60.0;
      m.mjd=-1.0;
      m.iodpoint=-1;
      redraw=1;
    }
    
    if (c=='l') {
      printf("Enter integration length (s): ");
      status=scanf("%f",&m.length);
      redraw=1;
    }

    // Toggle focal length
    if (c=='F') {
      fov++;
      if (fov>=sizeof(focallength)/sizeof(focallength[0]))
	fov=0;
      printf("Focallength: %.0f mm\n",focallength[fov]);
      m.fw=atan(0.5*6.3265/focallength[fov])*R2D;
      m.fh=atan(0.5*4.6389/focallength[fov])*R2D;
      //m.fw=atan(0.5*22.3/focallength[fov])*R2D;
      //m.fh=atan(0.5*14.9/focallength[fov])*R2D;
      printf("FOV: %.1fx%.1f\n",2*m.fw,2*m.fh);
      redraw=1;
    }

    if (c=='L') {
      if (Isatsel==0) {
	m.leoflag++;
	if (m.leoflag>2)
	  m.leoflag=-1;
	redraw=1;
      } else {
	printf("Unable, please reset satellite selection\n");
      }
    }
    
    // Find satellite
    if (c=='f') {
      printf("Enter NORAD Satellite number: ");
      status=scanf("%ld",&Isatsel);
      
      if (Isatsel!=0 && !print_tle(m.tlefile,Isatsel)) {
	printf("Satellite %ld not found!\n",Isatsel);
	Isatsel=-1;
      }      
      redraw=1;
    }
    
    // Measure
    if (c=='m') {
      if (strcmp(m.orientation,"horizontal")==0) {
	reverse(x,y,&azi,&alt);
	horizontal2equatorial(m.mjd,azi,alt,&ra,&de);
      } else if (strcmp(m.orientation,"equatorial")==0) {
	reverse(x,y,&ra,&de);
	equatorial2horizontal(m.mjd,ra,de,&azi,&alt);
      }
      
      printf("RA: %10.4f Dec: %10.4f Azi: %10.4f Alt: %10.4f\n%f %f\n",ra,de,modulo(azi-180.0,360.0),alt,x,y);
    }
    // Grid on/off
    if (c=='g' || c=='G') {
      if (m.grid==1)
	m.grid=0;
      else if (m.grid==0)
	m.grid=1;
      redraw=1;
    }
    
    // Exit
    if (c=='q' || c=='Q') {
      cpgend();
      exit(0);
    }
    
    // Recenter
    if (c=='c' || c=='C') {
      if (strcmp(m.orientation,"horizontal")==0) {
	reverse(x,y,&m.azi0,&m.alt0);
	horizontal2equatorial(m.mjd,m.azi0,m.alt0,&m.ra0,&m.de0);
      } else if (strcmp(m.orientation,"equatorial")==0) {
	reverse(x,y,&m.ra0,&m.de0);
	horizontal2equatorial(m.mjd,m.ra0,m.de0,&m.azi0,&m.alt0);
      }
      printf("Centered at: %8.4f %8.4f, %8.4f %8.4f\n",m.ra0,m.de0,m.azi0,m.alt0);
      redraw=1;
    }
    
    // Add to schedule
    if (c=='S') 
      schedule(m.nfd,m.ra0,m.de0);

    // Polar
    if (c=='z') {
      m.azi0=0.0;
      m.alt0=90.0;
      m.w=120.0;
      strcpy(m.orientation,"horizontal");
      m.level=1;
      redraw=1;
    }
    
    // South
    if (c=='s') {
      m.azi0=0.0;
      m.alt0=45.0;
      strcpy(m.orientation,"horizontal");
      m.level=3;
      redraw=1;
    }
    
    // North
    if (c=='n') {
      m.azi0=180.0;
      m.alt0=45.0;
      strcpy(m.orientation,"horizontal");
      m.level=3;
      redraw=1;
    }
    
    // East
    if (c=='e') {
      m.azi0=270.0;
      m.alt0=45.0;
      strcpy(m.orientation,"horizontal");
      m.level=3;
      redraw=1;
    }
    
    // West
    if (c=='w') {
      m.azi0=90.0;
      m.alt0=45.0;
      strcpy(m.orientation,"horizontal");
      m.level=3;
      redraw=1;
    }
    
    // Orientation
    if (c=='o' || c=='O') {
      if (strcmp(m.orientation,"horizontal")==0) {
	strcpy(m.orientation,"equatorial");
      } else if (strcmp(m.orientation,"equatorial")==0) {
	strcpy(m.orientation,"horizontal");
      }
      redraw=1;
    }
    
    // Level
    if (isdigit(c)) {
      m.level=c-'0';
      redraw=1;
    }
    
    // Zoom
    if (c=='-' && m.level>1) {
      m.level--;
      redraw=1;
    }
    if ((c=='+' || c=='=') && m.level<9) {
      m.level++;
      redraw=1;
    }
    
    // renew
    skymap_plot_renew();
  } 

  return 0;
}

// rotate vector
void rotate(int axis,float angle,float *x,float *y,float *z)
{
  float xx,yy,zz;
  
  if (axis==0) {
    xx= *x;
    yy= *y*cos(angle*D2R)- *z*sin(angle*D2R);
    zz= *z*cos(angle*D2R)+ *y*sin(angle*D2R);
  }
  if (axis==1) {
    xx= *x*cos(angle*D2R)- *z*sin(angle*D2R);
    yy= *y;
    zz= *z*cos(angle*D2R)+ *x*sin(angle*D2R);
  }
  if (axis==2) {
    xx= *x*cos(angle*D2R)- *y*sin(angle*D2R);
    yy= *y*cos(angle*D2R)+ *x*sin(angle*D2R);
    zz= *z;
  }
  
  *x=xx;
  *y=yy;
  *z=zz;

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
  x=3600.*fabs(x);
  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  hour=x;
  sec=floor(1000.0*sec)/1000.0;

  sprintf(date,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);

  return;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int f,int len)
{
  int i;
  double sec,deg,min;
  char sign;
  char *form[]={":: ",",, ","hms","   "};

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

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
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

// Convert Sexagesimal into Decimal
double sex2dec(char *s)
{
  double x;
  float deg,min,sec;
  char t[LIM];

  strcpy(t,s);

  deg=fabs(atof(strtok(t," :")));
  min=fabs(atof(strtok(NULL," :")));
  sec=fabs(atof(strtok(NULL," :")));

  x=(double) deg+(double) min/60.+(double) sec/3600.;
  if (s[0]=='-') x= -x;

  return x;
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

// Decode IOD Observations
struct observation decode_iod_observation(char *iod_line)
{
  int year,month,iday,hour,min;
  int format,epoch,me,xe,sign;
  int site_id;
  double sec,ra,mm,ss,de,dd,ds,day,mjd0;
  struct observation obs;
  char secbuf[6],sn[2],degbuf[3];

  // Strip newline
  iod_line[strlen(iod_line)-1]='\0';

  // Copy full line
  strcpy(obs.iod_line,iod_line);

  // Set usage
  obs.flag=1;

  // Get SSN
  sscanf(iod_line,"%5d",&obs.ssn);

  // Get site
  sscanf(iod_line+16,"%4d",&obs.site);

  // Decode date/time
  sscanf(iod_line+23,"%4d%2d%2d%2d%2d%5s",&year,&month,&iday,&hour,&min,secbuf);
  sec=atof(secbuf);
  sec/=pow(10,strlen(secbuf)-2);
  day=(double) iday+(double) hour/24.0+(double) min/1440.0+(double) sec/86400.0;
  obs.mjd=date2mjd(year,month,day);

  // Get uncertainty in time
  sscanf(iod_line+41,"%1d%1d",&me,&xe);
  obs.st=(float) me*pow(10,xe-8);

  // Skip empty observations
  if (strlen(iod_line)<64 || (iod_line[54]!='+' && iod_line[54]!='-')) 
    obs.flag=0;

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
  obs.sr=(float) me*pow(10,xe-8);
  
  // Decode position
  switch(format) 
    {
      // Format 1: RA/DEC = HHMMSSs+DDMMSS MX   (MX in seconds of arc)
    case 1 : 
      ra+=mm/60+ss/36000;
      de=sign*(de+dd/60+ds/3600);
      obs.sr/=3600.0;
      break;
      // Format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
    case 2:
      ra+=mm/60+ss/60000;
      de=sign*(de+dd/60+ds/6000);
      obs.sr/=60.0;
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
      printf("%s\n",iod_line);
      printf("IOD Format not implemented\n");
      obs.flag=0;
      break;
  }
  // Convert to degrees
  ra*=15.0;

  // Get precession epoch
  if (epoch==0) {
    obs.ra=ra;
    obs.de=de;
    return obs;
  } else if (epoch==4) {
    mjd0=33281.9235;
  } else if (epoch==5) {
    mjd0=51544.5;
  } else {
    printf("Observing epoch not implemented\n");
    obs.flag=0;
  }

  // Precess position
  precess(mjd0,ra,de,obs.mjd,&obs.ra,&obs.de);

  // Get horizontal position
  equatorial2horizontal(obs.mjd,obs.ra,obs.de,&obs.azi,&obs.alt);

  return obs;
}

void plot_iod(char *filename)
{
  int i=0;
  char line[LIM];
  FILE *file;
  struct observation obs;
  double azi,alt,rx,ry;
  float x,y;

  cpgsci(2);

  file=fopen(filename,"r");
  // Read data
  while (fgets(line,LIM,file)!=NULL) {
    if (strstr(line,"#")==NULL) {
      obs=decode_iod_observation(line);

      if (m.site_id==obs.site) {
	if (strcmp(m.orientation,"horizontal")==0) {
	  forward(obs.azi,obs.alt,&rx,&ry);
	} else if (strcmp(m.orientation,"equatorial")==0) {
	  forward(obs.ra,obs.de,&rx,&ry);
	}
	x=(float) rx;
	y=(float) ry;

	cpgpt1(x,y,4);
      }
    }
  }
  fclose(file);
  cpgsci(1);
  

  return;
}
