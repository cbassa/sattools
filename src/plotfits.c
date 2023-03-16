#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <wcslib/cel.h>
#include <cpgplot.h>
#include "qfits.h"
#include <gsl/gsl_multifit.h>
#include <getopt.h>

#define LIM 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define NMAX 4096

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
  float x0,y0;
  float a[3],b[3];
  double mjd;
} img;
struct catalog {
  int n;
  float x[NMAX],y[NMAX],mag[NMAX];
  double ra[NMAX],de[NMAX],rx[NMAX],ry[NMAX];
  int select[NMAX];
};
struct map {
  double lat,lng;
  float alt;
  int site_id;
  char observer[32];
} m;

struct image read_fits(char *filename,int pnum);
int fgetline(FILE *,char *,int);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double,double,double,double,double *,double *);
void lfit2d(float *x,float *y,float *z,int n,float *a);
struct catalog read_pixel_catalog(char *filename);
double gmst(double mjd);
double modulo(double x,double y);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);
double sex2dec(char *s);

// Read astrometric catalog
struct catalog read_astrometric_catalog(char *filename,float mmin,float sx,float sy,float angle)
{
  int i=0;
  FILE *file;
  char line[LIM];
  struct catalog c;
  double rx,ry,x,y,ra,de;
  struct star s;
  double d,dx,dy;
  double mjd0=51544.5;

  file=fopen(filename,"rb");
  if (file==NULL) {
    fprintf(stderr,"%s not found!\n",filename);
    exit(0);
  }
  while (!feof(file)) {
    fread(&s,sizeof(struct star),1,file);
    if (s.mag>mmin)
      continue;
    precess(mjd0,s.ra,s.de,img.mjd,&ra,&de);
    forward(img.ra0,img.de0,ra,de,&rx,&ry);
    x=img.x0+1.0/sx*(cos(angle*D2R)*rx+sin(angle*D2R)*ry);
    y=img.y0+1.0/sy*(-sin(angle*D2R)*rx+cos(angle*D2R)*ry);
      /*
    } else if (t.state==1) {
      dx=rx-t.a[0];
      dy=ry-t.b[0];
      d=t.a[1]*t.b[2]-t.a[2]*t.b[1];
      x=(t.b[2]*dx-t.a[2]*dy)/d;
      y=(t.a[1]*dy-t.b[1]*dx)/d;
    }
    */
    if (x>0.0 && x<img.naxis1 && y>0.0 && y<img.naxis2) {
      c.x[i]=x;
      c.y[i]=y;
      c.rx[i]=rx;
      c.ry[i]=ry;
      c.ra[i]=s.ra;
      c.de[i]=s.de;
      c.mag[i]=s.mag;
      c.select[i]=0;
      i++;
    }
  }
  fclose(file);
  c.n=i;

  return c;
}

// Read astrometric catalog
struct catalog reread_astrometric_catalog(char *filename,float mmin)
{
  int i=0;
  FILE *file;
  char line[LIM];
  struct catalog c;
  double rx,ry,x,y;
  struct star s;
  double d,dx,dy,ra,de;
  double mjd0=51544.5;

  file=fopen(filename,"rb");
  while (!feof(file)) {
    fread(&s,sizeof(struct star),1,file);
    if (s.mag>mmin)
      continue;
    precess(mjd0,s.ra,s.de,img.mjd,&ra,&de);
    forward(img.ra0,img.de0,ra,de,&rx,&ry);
    dx=rx-img.a[0];
    dy=ry-img.b[0];
    d=img.a[1]*img.b[2]-img.a[2]*img.b[1];
    x=(img.b[2]*dx-img.a[2]*dy)/d+img.x0;
    y=(img.a[1]*dy-img.b[1]*dx)/d+img.y0;
    if (x>0.0 && x<img.naxis1 && y>0.0 && y<img.naxis2) {
      c.x[i]=x;
      c.y[i]=y;
      c.rx[i]=rx;
      c.ry[i]=ry;
      c.ra[i]=s.ra;
      c.de[i]=s.de;
      c.mag[i]=s.mag;
      c.select[i]=0;
      i++;
    }
  }
  fclose(file);
  c.n=i;

  return c;
}

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

// Fit transformation
void fit_transformation(struct catalog cat,struct catalog ast,int nselect)
{
  int i,j;
  float *x,*y,*rx,*ry;

  x=(float *) malloc(sizeof(float)*nselect);
  y=(float *) malloc(sizeof(float)*nselect);
  rx=(float *) malloc(sizeof(float)*nselect);
  ry=(float *) malloc(sizeof(float)*nselect);
  
  for (i=0;i<nselect;i++) {
    for (j=0;j<cat.n;j++) {
      if (cat.select[j]==i+1) {
	x[i]=cat.x[j]-img.x0;
	y[i]=cat.y[j]-img.y0;
      }
    }
    for (j=0;j<ast.n;j++) {
      if (ast.select[j]==i+1) {
	rx[i]=ast.rx[j];
	ry[i]=ast.ry[j];
      } 
    }
  }
  
  lfit2d(x,y,rx,nselect,img.a);
  lfit2d(x,y,ry,nselect,img.b);

  

  return;
}

int match_catalogs(struct catalog *cat,struct catalog *ast,float rmax)
{
  int i,j,jmin,n,flag=0;
  float r,rmin;
  FILE *file;

  // Reset
  for (i=0;i<cat->n;i++)
    cat->select[i]=0;
  for (i=0;i<ast->n;i++)
    ast->select[i]=0;
  
  file=fopen("out.dat","w");
  for (i=0,n=0;i<cat->n;i++) {
    for (j=0,flag=0;j<ast->n;j++) {
      if (ast->select[j]!=0)
	continue;
      r=sqrt(pow(cat->x[i]-ast->x[j],2)+pow(cat->y[i]-ast->y[j],2));
      if (flag==0 || r<rmin) {
	rmin=r;
	jmin=j;
	flag=1;
      }
    }
    if (rmin<rmax) {
      fprintf(file,"%10.4f %10.4f %10.6f %10.6f\n",cat->x[i]-img.x0,cat->y[i]-img.y0,ast->ra[jmin],ast->de[jmin]);
      cat->select[i]=n+1;
      ast->select[jmin]=n+1;
      n++;
    }
  }
  fclose(file);

  printf("%d stars matched\n",n);
  return n;
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
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float x,y,r,rmin=1.0,rmax=10.0,mmin=5.0,mmax=10.0;
  struct catalog cat,ast;
  char c;
  int redraw=1,click=0,nselect=0,plotstars=1;
  char filename[128],sra[20],sde[20];
  float h,q,s=0.0,mag=9;
  FILE *file;
  char *env,starfile[128];

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(starfile,"%s/data/tycho2.dat",env);

  // Geographic position
  env=getenv("ST_COSPAR");
  get_site(atoi(env));


  // Read image
  img=read_fits(argv[1],0);
  sprintf(filename,"%s.cat",argv[1]);

  printf("Image read\n");

  // Initial transformation
  if (argc==7) {
    s=atof(argv[2]);
    img.ra0=atof(argv[3]);
    img.de0=atof(argv[4]);
    q=atof(argv[5]);
    mag=atof(argv[6]);
  } else {
    file=fopen("position.txt","r");
    if (file==NULL) {
      fprintf(stderr,"No position file found\n");
      return 0;
    }
    fscanf(file,"%s %s",sra,sde);
    fclose(file);
    
    // Get parameters
    img.ra0=15.0*sex2dec(sra);
    img.de0=sex2dec(sde);

    // Hour angle
    h=gmst(img.mjd)+m.lng-img.ra0;
    q=atan2(sin(h*D2R),(tan(m.lat*D2R)*cos(img.de0*D2R)-sin(img.de0*D2R)*cos(h*D2R)))*R2D;
    printf("Hour angle: %.3f deg, parallactic angle: %.3f deg\n",h,q);
  }
  img.x0=0.5*(float) img.naxis1;
  img.y0=0.5*(float) img.naxis2;

  // Read catalogs
  cat=read_pixel_catalog(filename);
  if (s==0.0)
    ast=read_astrometric_catalog(starfile,mag,-36.15,33.22,-q);
  else
    ast=read_astrometric_catalog(starfile,mag,-s,s,-q);

  // Plot image
  cpgopen("/xs");
  cpgwnad(0.0,img.naxis1,0.0,img.naxis2);
  cpgsfs(2);
  cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);

  // For ever loop
  for (;;) {      
    if (redraw==1) {
      cpgimag(img.z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,img.zmin,img.zmax,tr);
      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);
    
      // Plot catalogs
      if (plotstars==1) {
	cpgsci(3);
	for (i=0;i<cat.n;i++) {
	  if (cat.select[i]!=0) 
	    cpgpt1(cat.x[i],cat.y[i],6);
	  else
	    cpgpt1(cat.x[i],cat.y[i],4);
	}
      }
      cpgsci(4);
      for (i=0;i<ast.n;i++) {
	r=rmax-(rmax-rmin)*(ast.mag[i]-mmin)/(mmax-mmin);

	// Upscale for image size
	r*=img.naxis1/752.0;
	if (ast.select[i]!=0)
	  cpgpt1(ast.x[i],ast.y[i],6);
	cpgcirc(ast.x[i],ast.y[i],r);
      }
      cpgsci(1);
      redraw=0;
    }

    cpgcurs(&x,&y,&c);

    // Quit
    if (c=='q')
      break;

    // Fit
    if (c=='f' && nselect>=3) {
      fit_transformation(cat,ast,nselect);
      ast=reread_astrometric_catalog(starfile,mag+1);
      redraw=1;
    }

    // Reread
    if (c=='r') {
      ast=reread_astrometric_catalog(starfile,mag+1);
      redraw=1;
    }

    // Select pixel catalog
    if (c=='a' && click==0) {
      i=select_nearest(cat,x,y);
      cat.select[i]=nselect+1;
      redraw=1;
      click=1;
    }

    // Select catalog
    if (c=='b' && click==1) {
      i=select_nearest(ast,x,y);
      ast.select[i]=nselect+1;
      redraw=1;
      click=0;
      nselect++;
    }
    
    // 
    if (c=='p') {
      if (plotstars==1)
	plotstars=0;
      else if (plotstars==0)
	plotstars=1;
      redraw=1;
    }

    // Match catalogs
    if (c=='m') {
      nselect=match_catalogs(&cat,&ast,10.0);
      redraw=1;
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
  double s1,s2,avg,std;

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
  img.mjd=atof(qfits_query_hdr(filename,"MJD-OBS"));

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
  avg=s1/(float) (img.naxis1*img.naxis2);
  std=sqrt(s2/(float) (img.naxis1*img.naxis2)-avg*avg);
  printf("%f %f\n",avg,std);
  img.zmin=avg-4.0*std;
  img.zmax=avg+6.0*std;
  
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


// Linear 2D fit
void lfit2d(float *x,float *y,float *z,int n,float *a)
{
  int i;
  double chisq;
  gsl_matrix *X,*cov;
  gsl_vector *yy,*w,*c;

  X=gsl_matrix_alloc(n,3);
  yy=gsl_vector_alloc(n);
  w=gsl_vector_alloc(n);

  c=gsl_vector_alloc(3);
  cov=gsl_matrix_alloc(3,3);

  // Fill matrices
  for(i=0;i<n;i++) {
    gsl_matrix_set(X,i,0,1.0);
    gsl_matrix_set(X,i,1,x[i]);
    gsl_matrix_set(X,i,2,y[i]);
    
    gsl_vector_set(yy,i,z[i]);
    gsl_vector_set(w,i,1.0);
  }

  // Do fit
  gsl_multifit_linear_workspace *work=gsl_multifit_linear_alloc(n,3);
  gsl_multifit_wlinear(X,w,yy,c,cov,&chisq,work);
  gsl_multifit_linear_free(work);

  // Save parameters
  for (i=0;i<3;i++)
    a[i]=gsl_vector_get(c,(i));

  gsl_matrix_free(X);
  gsl_vector_free(yy);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return;
}

// Read pixel catalog
struct catalog read_pixel_catalog(char *filename)
{
  int i=0;
  FILE *file;
  char line[LIM];
  struct catalog c;

  // Read catalog
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"%s not found!\n",filename);
    exit(0);
  }
  while (fgetline(file,line,LIM)>0) {
    if (strstr(line,"#")!=NULL) 
      continue;
    sscanf(line,"%f %f %f",&c.x[i],&c.y[i],&c.mag[i]);
    c.select[i]=0;
    i++;
  }
  fclose(file);
  c.n=i;

  return c;
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
