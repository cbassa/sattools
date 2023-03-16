#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "cpgplot.h"
#include "sgdp4h.h"

#define LIM 80
#define NMAX 1024
#define MMAX 28368
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define FLAT (1.0/298.257)
#define XKMPAU 149597879.691 // AU in km

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

struct map {
  long satno;
  double lat,lng;
  double mjd;
  float alt,timezone;
  int length;
  char nfd[LIM],tlefile[LIM],observer[32];
  char datadir[LIM],tledir[LIM];
  int site_id;
  float l0,b0;
} m;
struct globe {
  int n;
  float l[MMAX],b[MMAX];
} glb;
struct sat {
  long Isat;
  double jd;
  double dx,dy,dz;
  double x,y,z,vx,vy,vz;
  double rsun,rearth;
  double psun,pearth,p;
  double r,ra,de;
  double azi,alt;
  double rx,ry;
};
void read_globe(void);
void plot_globe(void);
void initialize_setup(void);
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void mjd2date(double mjd,char *date,int length);
void usage();
void interactive_usage();
void nfd_now(char *s);
double gmst(double);
double dgmst(double);
double modulo(double,double);
void sunpos_xyz(double,xyz_t *,double *,double *);
void rotate(int axis,float angle,float *x,float *y,float *z);
void get_site(int site_id);

void plot_terminator(void)
{
  int i,j,j0,k,flag;
  xyz_t sunpos;
  double sra,sde,r,h;
  float l0,b0,l[NMAX+4],b[NMAX+4];
  float x,y,z;
  int isci;
  float theta,ang[]={0.0,-6.0,-12.0,-18.0};

  // Solar position
  sunpos_xyz(m.mjd,&sunpos,&sra,&sde);
  
  // GMST
  h=gmst(m.mjd);

  // Solar subpoint
  l0=modulo(sra-h,360.0);
  b0=sde;
  if (l0>180.0)
    l0-=360.0;
  
  // Loop over terminator boundaries
  for (k=0;k<4;k++) {
    for (i=0,j=0,flag=0;i<NMAX;i++,j++) {
      theta=2.0*M_PI*(float) i/(float) (NMAX-1);

      x=XKMPER*sin(ang[k]*D2R);
      y=XKMPER*sin(theta)*cos(ang[k]*D2R);
      z=XKMPER*cos(theta)*cos(ang[k]*D2R);

      rotate(1,b0,&x,&y,&z);
      rotate(2,l0,&x,&y,&z);
      
      r=sqrt(x*x+y*y+z*z);
      l[j]=atan2(y,x)*R2D;
      b[j]=asin(z/r)*R2D;
      l[j]=modulo(l[j],360.0);
      if (l[j]>180.0) 
	l[j]-=360.0;
      if (l[j]<-180.0) 
	l[j]+=360.0;
      
      // Passing limit left to right
      if (l[j]*l[j-1]<0.0 && fabs(l[j])>45.0 && flag==0 && k==0) {
	l[j+4]=l[j];
	b[j+4]=b[j];
	b[j]=b[j-1];
	b[j+3]=b[j-1];
	if (l[j-1]<l[j]) {
	  l[j]=-180.0;
	  l[j+1]=-180.0;
	  l[j+2]=180.0;
	  l[j+3]=180.0;
	} else {
	  l[j]=180.0;
	  l[j+1]=180.0;
	  l[j+2]=-180.0;
	  l[j+3]=-180.0;
	}
	if (b0<=0.0) {
	  b[j+1]=90.0;
	  b[j+2]=90.0;
	} else {
	  b[j+1]=-90.0;
	  b[j+2]=-90.0;
	}
	j+=4;
	flag=1;
      } 
    }
    
    if (k==0) {
      // Set night color
      cpgscr(16,0.0,0.0,0.2);

      // Plot night side
      cpgsci(16);
      cpgpoly(NMAX+4,l,b);

      // Plot terminator
      cpgsci(14);
      cpgline(NMAX+4,l,b);
      cpgsci(1);
    } else {
      // Plot twilight boundaries
      cpgsci(14);
      for (i=0,flag=0;i<NMAX;i++) {
	if (i>0 && l[i-1]*l[i]<0.0 && fabs(l[i-1]-l[i])>10.0)
	  flag=0;

	if (flag==0) {
	  cpgmove(l[i],b[i]);
	  flag=1;
	} else {
	  cpgdraw(l[i],b[i]);
	}
      }
      cpgsci(1);
    }
  }

  // Save sub solar position
  m.l0=l0;
  m.b0=b0;

  return;
}

void init_plot(char *psfile,float width,float aspect)
{
  
  cpgopen(psfile);
  cpgslw(2);
  //  cpgpap(width,aspect);
  cpgpap(0.0,aspect);

  return;
}

// Plot observing sites
void plot_sites(void)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64],text[8],filename[LIM];
  float isch;

  cpgqch(&isch);

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

    sprintf(text," %04d",id);
    cpgsci(2);
    cpgsch(0.5);
    cpgpt1(lng,lat,4);
    cpgtext(lng,lat,text);
    cpgsci(1);
  }
  fclose(file);
  cpgsch(isch);

  return;
}

// Plot observing sites
void plot_launch_sites(void)
{
  int i=0;
  char line[LIM];
  FILE *file;
  double lat,lng;
  char site[64],text[8],filename[LIM];
  float isch;

  cpgqch(&isch);

  sprintf(filename,"%s/data/launchsites.txt",m.datadir);
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
    sscanf(line,"%lf %lf",
	   &lat,&lng);
    strcpy(site,line+21);

    cpgsci(2);
    cpgsch(0.5);
    cpgpt1(lng,lat,4);
    cpgtext(lng,lat,site);
    cpgsci(1);
  }
  fclose(file);
  cpgsch(isch);

  return;
}


// Computes apparent position
struct sat apparent_position(double mjd)
{
  struct sat s;
  double jd,rsun,rearth;
  double dx,dy,dz;
  xyz_t satpos,obspos,satvel,sunpos;
  double sra,sde;

  // Sat ID
  s.Isat=Isat;

  // Get Julian Date
  jd=mjd+2400000.5;

  // Get positions
  satpos_xyz(jd,&satpos,&satvel);
  sunpos_xyz(mjd,&sunpos,&sra,&sde);

  // Sat positions
  s.x=satpos.x;
  s.y=satpos.y;
  s.z=satpos.z;
  s.vx=satvel.x;
  s.vy=satvel.y;
  s.vz=satvel.y;

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
  //  s.p=acos(((sunpos.x+satpos.x)*satpos.x+(sunpos.y+satpos.y)*satpos.y+(sunpos.z+satpos.z)*satpos.z)/(rsun*rearth))*R2D;

  s.p-=s.pearth;

  // Celestial position
  s.r=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
  s.ra=atan2(satpos.y,satpos.x)*R2D;
  s.de=asin(satpos.z/s.r)*R2D;
  
  return s;
}

// plot satellite track
void track_plot_track(char *tlefile,long satno,double mjd0)
{
  int i=0,nstep=500;
  orbit_t orb;
  xyz_t pos,vel;
  double jd,dt,h,l,b,l0,mjd;
  FILE *fp=NULL;
  float x,y,z,r,v;
  long imode;
  int isci;
  float isch;
  char norad[7];
  struct sat s;

  cpgqci(&isci);
  cpgqch(&isch);
  cpgsci(7);

  fp = fopen(tlefile, "rb");
  if(fp == NULL) {
    fatal_error("File open failed for reading \"%s\"", tlefile);
  }
  
  while(read_twoline(fp, satno, &orb) == 0) {
    //    print_orb(&orb);
    
    Isat = orb.satno;
    imode = init_sgdp4(&orb);
    
    if(imode == SGDP4_ERROR) continue;

    jd=mjd0+2400000.5;
    for (i=0;;i++) {
      //      if(satpos_xyz(jd, &pos, &vel) == SGDP4_ERROR) break;
      mjd=jd-2400000.5;
      s=apparent_position(mjd);

      h=gmst(mjd);

      x=s.x;
      y=s.y;
      z=s.z;
 
      // Celestial position
      r=sqrt(x*x+y*y+z*z);
      l=atan2(y,x)*R2D;
      b=asin(z/r)*R2D;
      l-=h;
      l=modulo(l,360.0);
      if (l>180.0) 
	l-=360.0;
      if (l<-180.0) 
	l+=360.0;
     
      // Visibility
      if (s.p<-s.psun)
	cpgsci(14);
      else if (s.p>-s.psun && s.p<s.psun)
	cpgsci(15);
      else if (s.p>s.psun)
	cpgsci(7);

      // Plot
      if (i==0) {
	sprintf(norad," %ld",Isat);
	cpgsch(0.6);
	cpgtext(l,b,norad);
	cpgsch(isch);
	cpgpt1(l,b,17);
	l0=l;
      }
      if (i==0 || fabs(l-l0)>10.0)
	cpgmove(l,b);
      else 
	cpgdraw(l,b);
      l0=l;

      // Do timestep
      r=sqrt(s.x*s.x+s.y*s.y+s.z*s.z);
      v=sqrt(s.vx*s.vx+s.vy*s.vy+s.vz*s.vz);
      dt=2.0*M_PI*r/(0.75*v*nstep);
      jd+=dt/86400.0;

      if (i==nstep)
	break;
    }
  }
  cpgsci(isci);
  cpgsch(isch);


  return;
}

void plot_map(void)
{
  int redraw=1;
  char text[256];
  float x,y;
  char c;

  for (;;) {
    if (redraw>0) {
      // Get present mjd
      if (m.mjd<0.0) {
	nfd_now(m.nfd);
	m.mjd=nfd2mjd(m.nfd);
      }

      cpgscr(0,0.0,0.0,0.0);
      cpgeras();

      // Create window
      cpgsvp(0.01,0.99,0.01,0.99);
      cpgwnad(-180.0,180.0,-90.0,90.0);

      // Set background
      cpgscr(0,0.0,0.0,0.5);
      cpgsci(0);
      cpgrect(-180.0,180.0,-90.0,90.0);
      cpgsci(1);
      cpgscr(0,0.0,0.0,0.0);
      cpgbox("BC",0.,0,"BC",0.,0);  
    
      // Top left string
      cpgsch(0.8);
      mjd2date(m.mjd,m.nfd,0);
      sprintf(text,"%s UTC",m.nfd);
      cpgmtxt("T",0.6,0.0,0.0,text);
      
      // Bottom string
      sprintf(text,"l: %d s",m.length);
      cpgmtxt("B",1.0,0.0,0.0,text);
      cpgsch(1.0);

      // Plot terminator
      plot_terminator();
      cpgsci(14);
      cpgbox("ABCG",30.,3,"ABCG",30.,3);  
      cpgsci(1);

      // Plot globe
      plot_globe();
      cpgsci(1);
      cpgbox("BCTS",30.,3,"BCTS",30.,3);

      // Plot sites
      //      plot_sites();

      // Plot launch sites
      plot_launch_sites();

      // Plot satellites
      track_plot_track(m.tlefile,m.satno,m.mjd);

      // Plot sub solar position
      cpgsci(7);
      cpgpt1(m.l0,m.b0,17);
      cpgsci(1);
    }
    
    // Reset redraw
    redraw=0;

    // Get cursor
    cpgcurs(&x,&y,&c);

    // Help
    if (c=='h') {
      interactive_usage();

      continue;
    }

    // Redraw
    if (c=='r') {
      m.mjd=-1.0;
      m.length=60;
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
    
    // Exit
    if (c=='q' || c=='Q') {
      cpgend();
      exit(0);
    }
  }

  return;
}



int main(int argc,char *argv[])
{
  int arg=0;

  // Initialize setup
  initialize_setup();

  // Decode options
  while ((arg=getopt(argc,argv,"t:c:i:s:l:h"))!=-1) {
    switch (arg) {
      
    case 't':
      strcpy(m.nfd,optarg);
      m.mjd=nfd2mjd(m.nfd);
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
  read_globe();

  // Initialize plot
  init_plot("/xs",0,0.75);
  
  plot_map();

  cpgend();



  return 0;
}

void read_globe(void)
{
  int i,status;
  FILE *file;
  char filename[LIM];

  sprintf(filename,"%s/data/globe.dat",m.datadir);
  file=fopen(filename,"r");

  for (i=0;i<MMAX;i++) {
    status=fscanf(file,"%f %f",&glb.b[i],&glb.l[i]);
  }
  fclose(file);
  glb.n=MMAX;

  return;
}

void plot_globe(void)
{
  int i,flag;

  for (i=0,flag=0;i<glb.n;i++) {
    if (glb.b[i]==9999.0) {
      flag=0;
      continue;
    }
    if (flag==0) {
      cpgmove(glb.l[i],glb.b[i]);
      flag=1;
    } else {
      cpgdraw(glb.l[i],glb.b[i]);
    }
  }

  return;
}

// Initialize setup
void initialize_setup(void)
{
  char *env;

  // Default parameters
  m.satno=0;
  m.timezone=0.0;
  m.length=60;
  nfd_now(m.nfd);
  m.mjd=nfd2mjd(m.nfd);

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

void usage()
{
  printf("usage: satmap -c TLEFILE [-t TIMESTAMP] [-s COSPARID] [-i SATNO]\n");
  printf("                [-l LENGTH] [-h]\n");
}

void interactive_usage()
{
  printf("Interactive help:");
  printf("r    Redraw\n");
  printf("\n");
  printf("<    Divide the integration length by a facor of 2\n");
  printf(">    Multiply the integration length by a facor of 2\n");
  printf("\n");
  printf(",    Increase time (+integration_length in seconds /(1 day))\n");
  printf(".    Roll back the time\n");
  printf("\n");
  printf("h    this interactive help\n");
  printf("q/Q  Exit\n");
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
  x=3600.*fabs(x);
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

// rotate vector
void rotate(int axis,float angle,float *x,float *y,float *z)
{
  float xx,yy,zz;
  float ca,sa;

  ca=cos(angle*D2R);
  sa=sin(angle*D2R);
  
  if (axis==0) {
    xx= *x;
    yy= *y*ca- *z*sa;
    zz= *z*ca+ *y*sa;
  }
  if (axis==1) {
    xx= *x*ca- *z*sa;
    yy= *y;
    zz= *z*ca+ *x*sa;
  }
  if (axis==2) {
    xx= *x*ca- *y*sa;
    yy= *y*ca+ *x*sa;
    zz= *z;
  }
  
  *x=xx;
  *y=yy;
  *z=zz;

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
  char abbrev[3],observer[64];
  char filename[LIM];

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
