#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qfits.h"
#include "cpgplot.h"
#include "cel.h"
#include <gsl/gsl_multifit.h>

#define LIM 256
#define NMAX 8192
#define D2R M_PI/180.0
#define R2D 180.0/M_PI

struct star {
  double ra,de;
  float pmra,pmde;
  float mag;
};
struct image {
  int naxis1,naxis2,naxis3;
  float *z;
  float zmin,zmax;
  double ra0,de0;
  float avg,std;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  float exptime;
  double mjd;
  char nfd[32];
  int cospar;
} img;
struct catalog {
  int n;
  float x[NMAX],y[NMAX],mag[NMAX];
  double ra[NMAX],de[NMAX],rx[NMAX],ry[NMAX];
  int select[NMAX];
};
struct observation {
  int satno,cospar;
  char desig[16],conditions,behavior;
  double mjd,ra,de;
  float terr,perr,tmid;
  char nfd[32],pos[32];
  int epoch,type;
  char iod_line[80];
  float x[3],y[3];
  int state;
};
struct image read_fits(char *filename,int pnum);
int fgetline(FILE *file,char *s,int lim);
int select_nearest(struct catalog c,float x,float y);

void plot_objects(char *filename)
{
  int i;
  FILE *file;
  float x0,y0,x1,y1,texp;
  int id;
  char line[LIM],catalog[128],dummy[128],text[8];

  file=fopen(filename,"r");
  while (fgetline(file,line,LIM)>0) {
    sscanf(line,"%s %f %f %f %f %f %d %s",dummy,&x0,&y0,&x1,&y1,&texp,&id,catalog);

    cpgsci(0);
    if (strstr(catalog,"classfd")!=NULL)
      cpgsci(4);

    cpgmove(x0,y0);
    cpgdraw(x1,y1);
    cpgpt1(x0,y0,4);
    sprintf(text," %05d",id);
    cpgtext(x0,y0,text);
  }
  fclose(file);
  cpgsci(1);

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

// Compute Date from Julian Day
void mjd2date(double mjd,char *date)
{
  double f,jd,dday;
  int z,alpha,a,b,c,d,e;
  double year,month,day,hour,min;
  double sec,x,fsec;

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
  fsec=1000.0*(sec-floor(sec));
  sprintf(date,"%04d%02d%02d%02d%02d%02.0f%03.0f",(int) year,(int) month,(int) day,(int) hour,(int) min,floor(sec),fsec);

  return;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int type)
{
  int i;
  double sec,deg,min,fmin;
  char sign;

  sign=(x<0 ? '-' : '+');
  x=60.*fabs(x);

  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;
  if (type==0)
    fmin=1000.0*(min-floor(min));
  else
    fmin=100.0*(min-floor(min));

  if (type==0)
    sprintf(s,"%02.0f%02.0f%03.0f",deg,floor(min),fmin);
  else
    sprintf(s,"%c%02.0f%02.0f%02.0f",sign,deg,floor(min),fmin);

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

// MJD to DOY
double mjd2doy(double mjd,int *yr)
{
  int year,month,k=2;
  int day;
  double doy;
  char nfd[32];
  
  mjd2date(mjd,nfd);

  sscanf(nfd,"%04d",&year);
  sscanf(nfd+4,"%02d",&month);
  sscanf(nfd+6,"%02d",&day);

  if (year%4==0 && year%400!=0)
    k=1;

  doy=floor(275.0*month/9.0)-k*floor((month+9.0)/12.0)+day-30;

  *yr=year;

  return doy;
}

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min;
  double mjd,dday;
  float sec;

  sscanf(date,"'%04d-%02d-%02dT%02d:%02d:%f'",&year,&month,&day,&hour,&min,&sec);

  dday=day+hour/24.0+min/1440.0+sec/86400.0;
  mjd=date2mjd(year,month,dday);

  return mjd;
}

void format_iod_line(struct observation *obs)
{
  int mt,xt,mp,xp;
  char string[10];
  
  // Time format
  sprintf(string,"%7.1e",obs->terr);
  mt=string[0]-'0';
  xt=atoi(string+4)+8;

  // Position format
  if (obs->type==2) {
    sprintf(string,"%7.1e",obs->perr);
    mp=string[0]-'0';
    xp=atoi(string+4)+8;
  } else {
    printf("Position format not implemented!\n");
  }

  sprintf(obs->iod_line,"%05d %c%c %-6s %04d %c %-17s %d%d %d%d %-14s %d%d %c",
	  obs->satno,
	  obs->desig[0],obs->desig[1],
	  obs->desig+2,
	  obs->cospar,
	  obs->conditions,
	  obs->nfd,
	  mt,xt,
	  obs->type,obs->epoch,
	  obs->pos,
	  mp,xp,
	  obs->behavior);

  return;
}

void find_designation(int satno0,char *desig0)
{
  FILE *file;
  int satno;
  char desig[16];
  char *env,filename[128];

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/desig.txt",env);

  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"Designation file not found!\n");
    exit(0);
  }
  while (!feof(file)) {
    fscanf(file,"%d %s",&satno,desig);
    if (satno==satno0) {
      strcpy(desig0,desig);
      break;
    }
  }
  fclose(file);

  return;
}

void write_observation(struct observation obs)
{
  FILE *file;

  file=fopen("observations.txt","a");
  fprintf(file,"%s\n",obs.iod_line);
  fclose(file);

  printf("Observation written\n");

  return;
}

// Reduce point
void reduce_point(struct observation *obs,struct image img,float tmid,float x,float y)
{
  int iframe,k;
  double ra,de,rx,ry;
  float dx,dy,dt;
  double mjd;
  char nfd[32],sra[15],sde[15];

  // Transform position
  dx=x-img.x0;
  dy=y-img.y0;
  rx=img.a[0]+img.a[1]*dx+img.a[2]*dy;
  ry=img.b[0]+img.b[1]*dx+img.b[2]*dy;
  reverse(img.ra0,img.de0,rx,ry,&ra,&de);

  dec2sex(ra/15.0,sra,0);
  dec2sex(de,sde,1);

  // Get time
  mjd=img.mjd+(double) tmid/86400.0;
  mjd2date(mjd,nfd);
  
  // Copy
  strcpy(obs->nfd,nfd);
  sprintf(obs->pos,"%s%s",sra,sde);

  return;
}

int main(int argc,char *argv[])
{
  int i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  int redraw=1,plotobj=1,click=0,nselect=0;
  float x,y,width;
  char c,pixcat[LIM],text[128];
  struct catalog cat,ast;
  float sx,sy,q;
  char *env,idfile[128];
  float r,rmin=1.0,rmax=10.0,mmin=1.0,mmax=8.0,mag=8.0;
  struct observation obs;
  double mjd,doy;
  int year;
  float frac=0.5;
  float fx=0.5,fy=0.333;

  // Environment variables
  env=getenv("ST_DATADIR");

  // Read image
  img=read_fits(argv[1],0);
  sprintf(idfile,"%s.id",argv[1]);

  // Set image aspect
  fx=0.5;
  fy=0.25*img.naxis1/img.naxis2;

  // Default observation
  obs.satno=99999;
  strcpy(obs.desig,"99999U");
  obs.cospar=atoi(env);
  obs.conditions='G';
  strcpy(obs.nfd,"YYYYMMDDHHMMSSsss");
  obs.terr=0.2;
  strcpy(obs.pos,"HHMMmmm+DDMMmm");
  strcpy(obs.iod_line,"");
  obs.perr=0.1;
  obs.epoch=5;
  obs.type=2;
  obs.behavior='S';
  obs.state=0;
  obs.cospar=img.cospar;

  // Get fake designation
  mjd=nfd2mjd(img.nfd);
  doy=mjd2doy(mjd,&year);
  sprintf(obs.desig,"%02d%03.0lfA",year-2000,doy+500);

  // Open server
  cpgopen("/xs");
  //  cpgpap(0.,1.25);
  cpgask(0);
  cpgsch(0.8);

  // Default limits
  width=(img.naxis1>img.naxis2) ? img.naxis1 : img.naxis2;
  xmin=0.0;
  xmax=img.naxis1;
  ymin=0.0;
  ymax=img.naxis2;
  zmin=img.zmin;
  zmax=img.zmax;

  // Forever loop
  for (;;) {
    if (redraw==1) {
      cpgeras();
      cpgsci(1);
      
      cpgsvp(0.1,0.9,0.1,0.85);
      cpgwnad(xmin,xmax,ymin,ymax);
      cpglab("x (pix)","y (pix)"," ");
      cpgsfs(2);
      cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);

      sprintf(text,"UT Date: %.23s  COSPAR ID: %04d",img.nfd+1,img.cospar);
      cpgmtxt("T",6.0,0.0,0.0,text);
      sprintf(text,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')",img.ra0,img.xrms,img.de0,img.yrms);
      cpgmtxt("T",4.8,0.0,0.0,text);
      sprintf(text,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d",img.naxis1*sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1])/3600.0,img.naxis2*sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2])/3600.0,sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1]),sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2]));
      cpgmtxt("T",3.6,0.0,0.0,text);

      cpgimag(img.z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

      // Plot fit
      if (obs.state==1) {
	cpgsci(4);
	cpgpt1(obs.x[0],obs.y[0],4);
	cpgmove(obs.x[1],obs.y[1]);
	cpgdraw(obs.x[2],obs.y[2]);
	cpgsci(1);
      } else if (obs.state==2) {
	cpgsci(4);
	cpgpt1(obs.x[0],obs.y[0],4);
	cpgsci(1);
      }

      if (plotobj==1)
	plot_objects(idfile);

      format_iod_line(&obs);
      cpgmtxt("T",1.0,0.5,0.5,obs.iod_line);

      redraw=0;
    }

    // Get cursor
    cpgcurs(&x,&y,&c);

    // Quit
    if (c=='q')
      break;

    // Change fraction
    if (c=='e') {
      if (frac>0.49 && frac<0.51)
	frac=1.0;
      else if (frac>0.51)
	frac=0.0;
      else if (frac<0.49)
	frac=0.5;
      printf("Fraction: %.1f\n",frac);
    }

    if (c=='p' || c=='X') {
      if (plotobj==1)
	plotobj=0;
      else
	plotobj=1;
      redraw=1;
    }

    // Get designation
    if (c=='d') {
      printf("Provide satellite number: ");
      scanf("%d",&obs.satno);
      find_designation(obs.satno,obs.desig);
      redraw=1;
      continue;
    }
    // Center
    if (c=='c') {
      xmin=x-fx*width;
      xmax=x+fx*width;
      ymin=y-fy*width;
      ymax=y+fy*width;
      redraw=1;
      continue;
    }

    if (isdigit(c)) {
      width=1000.0/(c-'0');
      xmin=x-fx*width;
      xmax=x+fx*width;
      ymin=y-fy*width;
      ymax=y+fy*width;
      redraw=1;
      continue;
    }



    // Zoom
    if (c=='z' || c=='+') {
      width/=1.5;
      xmin=x-fx*width;
      xmax=x+fx*width;
      ymin=y-fy*width;
      ymax=y+fy*width;
      redraw=1;
      continue;
    }

    // Unzoom
    if (c=='x' || c=='+' || c=='=') {
      width*=1.5;
      xmin=x-fx*width;
      xmax=x+fx*width;
      ymin=y-fy*width;
      ymax=y+fy*width;
      redraw=1;
      continue;
    }

    // Reset
    if (c=='r') {
      width=(img.naxis1>img.naxis2) ? img.naxis1 : img.naxis2;
      xmin=0.0;
      xmax=img.naxis1;
      ymin=0.0;
      ymax=img.naxis2;
      redraw=1;
      continue;
    }

    // Write obs
    if (c=='w') {
      write_observation(obs);
      continue;
    }

    // Measure
    if (c=='m' || c=='D') {
      reduce_point(&obs,img,frac*img.exptime,x,y);
      obs.x[0]=x;
      obs.y[0]=y;
      obs.state=2;
      redraw=1;
      continue;
    }
  }

  cpgend();


  return 0;
}

// Read fits image
struct image read_fits(char *filename,int pnum)
{
  int i,j,k,l,m;
  qfitsloader ql;
  char key[FITS_LINESZ+1] ;
  struct image img;
  float s1,s2,avg,std;

  // Set plane
  ql.xtnum = 0;
  ql.pnum = pnum;

  // Set loadtype
  ql.ptype = PTYPE_FLOAT;

  // Set filename
  ql.filename=filename;

  // Image size
  img.naxis1=atoi(qfits_query_hdr(filename,"NAXIS1"));
  img.naxis2=atoi(qfits_query_hdr(filename,"NAXIS2"));

  // MJD
  img.mjd=atof(qfits_query_hdr(filename,"MJD-OBS"));
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

  // Initialize load
  if (qfitsloader_init(&ql) != 0) 
    printf("Error initializing data loading\n");

  // Test load
  if (qfits_loadpix(&ql) != 0) 
    printf("Error loading actual data\n");

  // Allocate image memory
  img.z=(float *) malloc(sizeof(float) * img.naxis1*img.naxis2);

  // Fill z array
  for (i=0,l=0,m=0;i<img.naxis1;i++) {
    for (j=0;j<img.naxis2;j++) {
      img.z[l]=ql.fbuf[l];
      l++;
    }
  }

  // Get levels
  for (i=0,s1=0.0,s2=0.0;i<img.naxis1*img.naxis2;i++) {
    s1+=img.z[i];
    s2+=img.z[i]*img.z[i];
  }
  img.avg=s1/(float) (img.naxis1*img.naxis2);
  img.std=sqrt(s2/(float) (img.naxis1*img.naxis2)-img.avg*img.avg);
  img.zmin=img.avg-4.0*img.std;
  img.zmax=img.avg+12.0*img.std;
  
  return img;
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


// Select nearest object
int select_nearest(struct catalog c,float x,float y)
{
  int i,imin;
  float r,rmin;

  for (i=0;i<c.n;i++) {
    r=sqrt(pow(x-c.x[i],2)+pow(y-c.y[i],2));
    if (i==0 || r<rmin) {
      imin=i;
      rmin=r;
    }
  }

  return imin;
}

