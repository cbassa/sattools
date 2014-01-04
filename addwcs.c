#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cel.h"
#include "cpgplot.h"
#include "qfits.h"
#include <gsl/gsl_multifit.h>
#include <getopt.h>

#define NMAX 8192
#define LIM 128
#define D2R M_PI/180.0
#define R2D 180.0/M_PI

struct image {
  int naxis,naxis1,naxis2,nframes;
  float *zavg,*zstd,*zmax,*znum;
  double ra0,de0;
  float x0,y0;
  float a[2],b[2];
  double mjd;
  float *dt;
};
struct transformation {
  double mjd;
  double ra0,de0;
  float a[3],b[3];
  float x0,y0;
  float xrms,yrms,rms;
};
struct star {
  double ra,de;
  float pmra,pmde;
  float mag;
};
struct catalog {
  int n;
  float x[NMAX],y[NMAX],imag[NMAX],fm[NMAX],fb[NMAX],bg[NMAX];
  double ra[NMAX],de[NMAX],vmag[NMAX];
  double rx[NMAX],ry[NMAX];
  float xres[NMAX],yres[NMAX],res[NMAX];
  int usage[NMAX];
};
struct image read_fits(char *filename);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double ra0,double de0,double x,double y,double *ra,double *de);
double gmst(double mjd);
double modulo(double x,double y);
int fgetline(FILE *file,char *s,int lim);
struct catalog match_catalogs(char *pixcat,char *astcat,struct transformation t,struct image img,float rmax,float mmin);
void plot_astrometric_catalog(struct transformation t,struct image img,float mmin);
void plot_pixel_catalog(char *filename);
void lfit2d(float *x,float *y,float *z,int n,float *a);
void add_fits_keywords(struct transformation t,char *filename);
void modify_fits_keywords(struct transformation t,char *filename);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);

void plot_image(struct image img,struct transformation t,struct catalog c,char *filename,float mmin)
{
  int i;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[]={0.0,0.2,0.4,0.6,1.0};
  float heat_r[]={0.0,0.5,1.0,1.0,1.0};
  float heat_g[]={0.0,0.0,0.5,1.0,1.0};
  float heat_b[]={0.0,0.0,0.0,0.3,1.0};
  float zmin,zmax,zavg,zstd;

  for (i=0,zavg=0.0;i<img.naxis1*img.naxis2;i++)
    zavg+=img.zavg[i];
  zavg/=(float) img.naxis1*img.naxis2;
  for (i=0,zstd=0.0;i<img.naxis1*img.naxis2;i++)
    zstd+=pow(img.zavg[i]-zavg,2);
  zstd=sqrt(zstd/(float) (img.naxis1*img.naxis2));
  zmin=zavg-2*zstd;
  zmax=zavg+6*zstd;

  cpgopen("1/xs");
  cpgwnad(0.0,img.naxis1,0.0,img.naxis2);
  cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
    
  cpgimag(img.zavg,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
  cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

  cpgsci(3);
  plot_pixel_catalog(filename);

  cpgsci(4);
  plot_astrometric_catalog(t,img,mmin);
  cpgsci(2);
  for (i=0;i<c.n;i++)
    cpgpt1(c.x[i]+t.x0,c.y[i]+t.y0,24);

  cpgend();

  return;
}

void usage(float mmin,float rmin)
{
  printf("addwcs: Add/fit World Coordinate System to a FITS file\n\n");

  printf("-f <file>:    FITS file to add/fit WCS to [required]\n");
  printf("-r <file>:    FITS file with reference WCS [required]\n");
  printf("-m <float>:   Magnitude cut-off in Tycho-2 catalog [optional; default %.1f]\n",mmin);
  printf("-R <float>:   Radius cut-off for matching [optional; default %.1f pix]\n",rmin);
  printf("-p            Plot image and selected stars [optional]\n");
  printf("-a            Add WCS keywords to input file (instead of modify) [optional]\n");
  printf("-t            Track on a fixed RA/Dec (correct for field rotation)\n");
  printf("-h            Print this help\n");

  return;
}

// Get reference transformation
struct transformation reference(char *filename)
{
  struct transformation t;

  t.mjd=atof(qfits_query_hdr(filename,"MJD-OBS"));
  t.ra0=atof(qfits_query_hdr(filename,"CRVAL1"));
  t.de0=atof(qfits_query_hdr(filename,"CRVAL2"));
  t.x0=atof(qfits_query_hdr(filename,"CRPIX1"));
  t.y0=atof(qfits_query_hdr(filename,"CRPIX2"));
  t.a[0]=0.0;
  t.a[1]=3600.0*atof(qfits_query_hdr(filename,"CD1_1"));
  t.a[2]=3600.0*atof(qfits_query_hdr(filename,"CD1_2"));
  t.b[0]=0.0;
  t.b[1]=3600.0*atof(qfits_query_hdr(filename,"CD2_1"));
  t.b[2]=3600.0*atof(qfits_query_hdr(filename,"CD2_2"));

  return t;
}

void rotate(float theta,float *x,float *y)
{
  float ct,st;
  float x0,y0;
  
  ct=cos(theta*D2R);
  st=sin(theta*D2R);
  x0= *x;
  y0= *y;

  *x=ct*x0-st*y0;
  *y=st*x0+ct*y0;

  return;
}

int main(int argc,char *argv[])
{
  int i,j,k,l,m;
  struct transformation t;
  struct image img;
  char *fitsfile=NULL,*reffile=NULL,catfile[128],calfile[128];
  FILE *outfile;
  struct catalog c;
  float mmin=10.0,rmin=10.0;
  double mjd0=51544.5,ra0,de0,ra1,de1;
  float q0,q1;
  float rmsmin;
  float x[NMAX],y[NMAX],rx[NMAX],ry[NMAX];
  int arg=0,plot=0,add=0,track=0;
  char *env,starfile[128];

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(starfile,"%s/data/tycho2.dat",env);

  // Decode options
  if (argc>1) {
    while ((arg=getopt(argc,argv,"f:r:m:R:hpnta"))!=-1) {
      switch (arg) {
	
      case 'f':
	fitsfile=optarg;
	break;
	
      case 'r':
	reffile=optarg;
	break;
	
      case 'm':
	mmin=atof(optarg);
	break;

      case 't':
	track=1;
	break;

      case 'R':
	rmin=atof(optarg);
	break;
	
      case 'p':
	plot=1;
	break;

      case 'a':
	add=1;
	break;
	
      case 'h':
	usage(mmin,rmin);
	return 0;
	
      default:
	usage(mmin,rmin);
	return 0;
      }
    } 
  } else {
    usage(mmin,rmin);
    return 0;
  }

  // Check if minimum input is provided
  if (fitsfile==NULL || reffile==NULL) {
    usage(mmin,rmin);
    return 0;
  }

  // Check this is indeed a FITS file 
  if (is_fits_file(fitsfile)!=1) {
    printf("%s is not a FITS file\n",fitsfile);
    return -1 ;
  }
    
  // Check this is indeed a FITS file 
  if (is_fits_file(reffile)!=1) {
    printf("%s is not a FITS file\n",reffile);
    return -1 ;
  }  

  // Read fits file
  img=read_fits(fitsfile);
  sprintf(catfile,"%s.cat",fitsfile);
  sprintf(calfile,"%s.cal",fitsfile);

  // Read reference transformation
  t=reference(reffile);

  // Correct astrometry for fixed or tracked setup
  if (track==0) {
    precess(mjd0,t.ra0,t.de0,t.mjd,&ra1,&de1);
    ra1=modulo(ra1+gmst(img.mjd)-gmst(t.mjd),360.0);
    precess(img.mjd,ra1,de1,mjd0,&t.ra0,&t.de0);
  }

  // Match catalog
  c=match_catalogs(catfile,starfile,t,img,rmin,mmin);

  
  // Plot
  if (plot==1)
    plot_image(img,t,c,catfile,mmin);

  // Do fit
  if (c.n>10) {
    for (l=0;l<10;l++) {
      for (j=0;j<5;j++) {
	// Transform
	for (i=0;i<c.n;i++) 
	  forward(t.ra0,t.de0,c.ra[i],c.de[i],&c.rx[i],&c.ry[i]);
      
	// Select
	for (i=0,k=0;i<c.n;i++) {
	  if (c.usage[i]==1) {
	    x[k]=c.x[i];
	    y[k]=c.y[i];
	    rx[k]=(float) c.rx[i];
	    ry[k]=(float) c.ry[i];
	    k++;
	  }
	}
	// Fit
	lfit2d(x,y,rx,k,t.a);
	lfit2d(x,y,ry,k,t.b);
	
	// Move reference point
	reverse(t.ra0,t.de0,t.a[0],t.b[0],&ra0,&de0);
	t.ra0=ra0;
	t.de0=de0;
      }
      
      // Compute and plot residuals
      for (i=0,t.xrms=0.0,t.yrms=0.0,m=0;i<c.n;i++) {
	if (c.usage[i]==1) {
	  c.xres[i]=c.rx[i]-(t.a[0]+t.a[1]*c.x[i]+t.a[2]*c.y[i]);
	  c.yres[i]=c.ry[i]-(t.b[0]+t.b[1]*c.x[i]+t.b[2]*c.y[i]);
	  
	  c.res[i]=sqrt(c.xres[i]*c.xres[i]+c.yres[i]*c.yres[i]);
	  t.xrms+=c.xres[i]*c.xres[i];
	  t.yrms+=c.yres[i]*c.yres[i];
	  t.rms+=c.xres[i]*c.xres[i]+c.yres[i]*c.yres[i];
	  m++;
	}
      }
      t.xrms=sqrt(t.xrms/(float) m);
      t.yrms=sqrt(t.yrms/(float) m);
      t.rms=sqrt(t.rms/(float) m);
      
      // Deselect outliers
      for (i=0;i<c.n;i++) {
	if (c.res[i]>2*t.rms)
	  c.usage[i]=0;
      }
    }
  } else {
    t.xrms=0.0;
    t.yrms=0.0;
    t.rms=0.0;
  }

  // Print results
  outfile=fopen(calfile,"w");
  for (i=0;i<c.n;i++) 
    if (c.usage[i]==1)
      fprintf(outfile,"%10.4f %10.4f %10.6f %10.6f %8.3f %8.3f %8.3f %8.3f %8.3f\n",c.x[i],c.y[i],c.ra[i],c.de[i],c.vmag[i],c.imag[i],c.fb[i],c.fm[i],c.bg[i]);
  fclose(outfile);

  printf("%s %8.4lf %8.4lf ",fitsfile,t.ra0,t.de0);
  printf("%3d/%3d %6.1f %6.1f %6.1f\n",m,c.n,t.xrms,t.yrms,t.rms);

  // Add keywords
  if (add==1)
    add_fits_keywords(t,fitsfile);
  else
    modify_fits_keywords(t,fitsfile);

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

  // Image size
  img.naxis=atoi(qfits_query_hdr(filename,"NAXIS"));
  img.naxis1=atoi(qfits_query_hdr(filename,"NAXIS1"));
  img.naxis2=atoi(qfits_query_hdr(filename,"NAXIS2"));

  // MJD
  img.mjd=(double) atof(qfits_query_hdr(filename,"MJD-OBS"));

  // Set parameters
  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;

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
      ql.pnum=k;;
      
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

// Match catalogs
struct catalog match_catalogs(char *pixcat,char *astcat,struct transformation t,struct image img,float rmax,float mmin)
{
  int i=0,imin,j,k,np;
  FILE *file;
  char line[LIM];
  struct star s;
  double rx,ry,d,dx,dy;
  int usage[NMAX];
  float xp[NMAX],yp[NMAX],mp[NMAX],x,y,fb[NMAX],fm[NMAX],bg[NMAX];
  struct catalog c;
  float r,rmin;

  // Read pixel catalog
  file=fopen(pixcat,"r");
  if (file==NULL) {
    printf("pixel catalog not found\n");
    exit(-1);
  }
  while (fgetline(file,line,LIM)>0) {
    if (strstr(line,"#")!=NULL) 
      continue;
    sscanf(line,"%f %f %f %f %f %f",&xp[i],&yp[i],&mp[i],&fb[i],&fm[i],&bg[i]);
    usage[i]=1;
    i++;
  }
  fclose(file);
  np=i;

  // Denominator
  d=t.a[1]*t.b[2]-t.a[2]*t.b[1];

  // Read astrometric catalog
  file=fopen(astcat,"rb");
  if (file==NULL) {
    printf("astrometric catalog not found\n");
    exit(-1);
  }
  j=0;
  while (!feof(file)) {
    fread(&s,sizeof(struct star),1,file);
    if (s.mag>mmin)
      continue;
    r=acos(sin(t.de0*D2R)*sin(s.de*D2R)+cos(t.de0*D2R)*cos(s.de*D2R)*cos((t.ra0-s.ra)*D2R))*R2D;
    if (r>90.0)
      continue;
    forward(t.ra0,t.de0,s.ra,s.de,&rx,&ry);

    dx=rx-t.a[0];
    dy=ry-t.b[0];
    x=(t.b[2]*dx-t.a[2]*dy)/d+t.x0;
    y=(t.a[1]*dy-t.b[1]*dx)/d+t.y0;

    // On image
    if (x>0.0 && x<img.naxis1 && y>0.0 && y<img.naxis2) {
      // Loop over pixel catalog
      for (i=0;i<np;i++) {
	r=sqrt(pow(xp[i]-x,2)+pow(yp[i]-y,2));
	if (i==0 || r<rmin) {
	  rmin=r;
	  imin=i;
	}
      }

      // Select
      if (rmin<rmax && usage[imin]==1) {
	c.x[j]=xp[imin]-t.x0;
	c.y[j]=yp[imin]-t.y0;
	c.imag[j]=mp[imin];
	c.fb[j]=fb[imin];
	c.fm[j]=fm[imin];
	c.bg[j]=bg[imin];
	c.ra[j]=s.ra;
	c.de[j]=s.de;
	c.vmag[j]=s.mag;
	c.usage[j]=1;
	usage[imin]=0;
	j++;
      }
    }
  }
  fclose(file);
  c.n=j;

  return c;
}

// Plot astrometric catalog
void plot_astrometric_catalog(struct transformation t,struct image img,float mmin)
{
  int i=0;
  FILE *file;
  struct star s;
  double rx,ry,d,r;
  double ra,de;
  float x,y;
  char *env,starfile[128];

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(starfile,"%s/data/tycho2.dat",env);

  d=t.a[1]*t.b[2]-t.a[2]*t.b[1];

  file=fopen(starfile,"rb");
  while (!feof(file)) {
    fread(&s,sizeof(struct star),1,file);
    if (s.mag>mmin)
      continue;
    r=acos(sin(t.de0*D2R)*sin(s.de*D2R)+cos(t.de0*D2R)*cos(s.de*D2R)*cos((t.ra0-s.ra)*D2R))*R2D;
    if (r>90.0)
      continue;
    forward(t.ra0,t.de0,s.ra,s.de,&rx,&ry);
    x=(t.b[2]*rx-t.a[2]*ry)/d+t.x0;
    y=(t.a[1]*ry-t.b[1]*rx)/d+t.y0;
    if (x>0.0 && x<img.naxis1 && y>0.0 && y<img.naxis2)
      cpgpt1(x,y,24);
  }
  fclose(file);

  return;
}

// Plot pixel catalog
void plot_pixel_catalog(char *filename)
{
  int i=0;
  FILE *file;
  char line[LIM];
  float x,y,mag;

  // Read catalog
  file=fopen(filename,"r");
  while (fgetline(file,line,LIM)>0) {
    if (strstr(line,"#")!=NULL) 
      continue;
    sscanf(line,"%f %f %f",&x,&y,&mag);

    cpgpt1(x,y,4);

    i++;
  }
  fclose(file);

  return;
}

// Linear 2D fit
void lfit2d(float *x,float *y,float *z,int n,float *a)
{
  int i,j,m;
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

// Add FITS keywords
void add_fits_keywords(struct transformation t,char *filename) 
{
  int i,j,k,l,m;
  int naxis1,naxis2,naxis3;
  qfits_header *qh;
  qfitsdumper qd;
  qfitsloader ql;
  char key[FITS_LINESZ+1];
  char val[FITS_LINESZ+1];
  char com[FITS_LINESZ+1];
  char lin[FITS_LINESZ+1];
  FILE *file;
  float *fbuf;

  naxis1=atoi(qfits_query_hdr(filename,"NAXIS1"));
  naxis2=atoi(qfits_query_hdr(filename,"NAXIS2"));
  naxis3=atoi(qfits_query_hdr(filename,"NAXIS3"));

  fbuf=malloc(sizeof(float)*naxis1*naxis2*naxis3);

  // Read header
  qh=qfits_header_read(filename);

  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;
  for (k=0,l=0;k<naxis3;k++) {
    ql.pnum=k;
    // Initialize load
    if (qfitsloader_init(&ql) != 0) 
      printf("Error initializing data loading\n");

    // Test load
    if (qfits_loadpix(&ql) != 0) 
      printf("Error loading actual data\n");

    for (i=0,m=0;i<naxis1;i++) {
      for (j=0;j<naxis2;j++) {
	fbuf[l]=ql.fbuf[m];
	l++;
	m++;
      }
    }
  }

  sprintf(val,"%e",t.yrms/3600.0);
  qfits_header_add_after(qh,"MJD-OBS","CRRES2",val," ",NULL);
  sprintf(val,"%e",t.xrms/3600.0);
  qfits_header_add_after(qh,"MJD-OBS","CRRES1",val," ",NULL);
  qfits_header_add_after(qh,"MJD-OBS","CUNIT2","'deg'"," ",NULL);
  qfits_header_add_after(qh,"MJD-OBS","CUNIT1","'deg'"," ",NULL);
  qfits_header_add_after(qh,"MJD-OBS","CTYPE2","'DEC--TAN'"," ",NULL);
  qfits_header_add_after(qh,"MJD-OBS","CTYPE1","'RA---TAN'"," ",NULL);
  sprintf(val,"%e",t.b[2]/3600.0);
  qfits_header_add_after(qh,"MJD-OBS","CD2_2",val," ",NULL);
  sprintf(val,"%e",t.b[1]/3600.0);
  qfits_header_add_after(qh,"MJD-OBS","CD2_1",val," ",NULL);
  sprintf(val,"%e",t.a[2]/3600.0);
  qfits_header_add_after(qh,"MJD-OBS","CD1_2",val," ",NULL);
  sprintf(val,"%e",t.a[1]/3600.0);
  qfits_header_add_after(qh,"MJD-OBS","CD1_1",val," ",NULL);
  sprintf(val,"%f",t.de0);
  qfits_header_add_after(qh,"MJD-OBS","CRVAL2",val," ",NULL);
  sprintf(val,"%f",t.ra0);
  qfits_header_add_after(qh,"MJD-OBS","CRVAL1",val," ",NULL);
  sprintf(val,"%f",t.y0);
  qfits_header_add_after(qh,"MJD-OBS","CRPIX2",val," ",NULL);
  sprintf(val,"%f",t.x0);
  qfits_header_add_after(qh,"MJD-OBS","CRPIX1",val," ",NULL);

  file=fopen(filename,"w");
  qfits_header_dump(qh,file);
  fclose(file);

  qfits_header_destroy(qh);

  qd.filename=filename;
  qd.npix=naxis1*naxis2*naxis3;
  qd.ptype=PTYPE_FLOAT;
  qd.fbuf=fbuf;
  qd.out_ptype=-32;

  qfits_pixdump(&qd);
  free(fbuf);

  return;
}

// Modify FITS keywords
void modify_fits_keywords(struct transformation t,char *filename) 
{
  char card[FITS_LINESZ+1];
  char key[FITS_LINESZ+1];
  char val[FITS_LINESZ+1];
  char com[FITS_LINESZ+1];

  sprintf(val,"%f",t.x0);
  keytuple2str(card,"CRPIX1",val,"");
  qfits_replace_card(filename,"CRPIX1",card);

  sprintf(val,"%f",t.y0);
  keytuple2str(card,"CRPIX2",val,"");
  qfits_replace_card(filename,"CRPIX2",card);

  sprintf(val,"%f",t.ra0);
  keytuple2str(card,"CRVAL1",val,"");
  qfits_replace_card(filename,"CRVAL1",card);

  sprintf(val,"%f",t.de0);
  keytuple2str(card,"CRVAL2",val,"");
  qfits_replace_card(filename,"CRVAL2",card);

  sprintf(val,"%e",t.a[1]/3600.0);
  keytuple2str(card,"CD1_1",val,"");
  qfits_replace_card(filename,"CD1_1",card);

  sprintf(val,"%e",t.a[2]/3600.0);
  keytuple2str(card,"CD1_2",val,"");
  qfits_replace_card(filename,"CD1_2",card);

  sprintf(val,"%e",t.b[1]/3600.0);
  keytuple2str(card,"CD2_1",val,"");
  qfits_replace_card(filename,"CD2_1",card);

  sprintf(val,"%e",t.b[2]/3600.0);
  keytuple2str(card,"CD2_2",val,"");
  qfits_replace_card(filename,"CD2_2",card);

  sprintf(val,"%e",t.xrms/3600.0);
  keytuple2str(card,"CRRES1",val,"");
  qfits_replace_card(filename,"CRRES1",card);

  sprintf(val,"%e",t.yrms/3600.0);
  keytuple2str(card,"CRRES2",val,"");
  qfits_replace_card(filename,"CRRES2",card);


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
