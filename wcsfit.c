#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cel.h"
#include "cpgplot.h"
#include "qfits.h"
#include <gsl/gsl_multifit.h>

#define LIM 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define NMAX 2048

struct catalog {
  int n;
  float x[NMAX],y[NMAX];
  double ra[NMAX],de[NMAX];
  float rx[NMAX],ry[NMAX];
  float xres[NMAX],yres[NMAX],res[NMAX];
  float xrms,yrms,rms;
  int usage[NMAX];
};
struct image {
  int naxis1,naxis2,nframes;
  float *zavg,*zstd,*zmax,*znum;
  double ra0,de0;
  float x0,y0;
  float a[2],b[2];
  double mjd;
  float *dt;
};
struct transformation {
  double ra0,de0;
  float a[3],b[3];
  float x0,y0;
};
int fgetline(FILE *file,char *s,int lim);
void forward(double ra0,double de0,double ra,double de,float *x,float *y);
void reverse(double ra0,double de0,float x,float y,double *ra,double *de);
struct catalog read_catalog(char *filename);
void lfit2d(float *x,float *y,float *z,int n,float *a);
void add_fits_keywords(struct transformation t,char *filename);
struct image read_fits(char *filename);

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

  return;
}


int main(int argc,char *argv[])
{
  int i,j,k,l,m;
  struct catalog c;
  struct transformation t;
  double ra0,de0;
  float rmsmin;
  float x[NMAX],y[NMAX],rx[NMAX],ry[NMAX];
  struct image img;
  char filename[128];

  if (argc==1) 
    strcpy(filename,"test.fits");
  else if (argc==2)
    strcpy(filename,argv[1]);

  img=read_fits(filename);
  printf("files read\n");
  c=read_catalog("out.dat");
  printf("files read\n");

  // Initial fit
  t.ra0=c.ra[0];
  t.de0=c.de[0];
  t.x0=(float) img.naxis1/2.0;
  t.y0=(float) img.naxis2/2.0;

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
	  rx[k]=c.rx[i];
	  ry[k]=c.ry[i];
	  k++;
	}
      }

      // Fit
      lfit2d(x,y,rx,k,t.a);
      lfit2d(x,y,ry,k,t.b);
      printf("%f %f %f %f %f %f %f %f\n",t.ra0,t.de0,t.a[0],t.a[1],t.a[2],t.b[0],t.b[1],t.b[2]);
      
      // Move reference point
      reverse(t.ra0,t.de0,t.a[0],t.b[0],&ra0,&de0);
      t.ra0=ra0;
      t.de0=de0;
    }

    // Compute and plot residuals
    for (i=0,c.xrms=0.0,c.yrms=0.0,m=0;i<c.n;i++) {
      if (c.usage[i]==1) {
	c.xres[i]=c.rx[i]-(t.a[0]+t.a[1]*c.x[i]+t.a[2]*c.y[i]);
	c.yres[i]=c.ry[i]-(t.b[0]+t.b[1]*c.x[i]+t.b[2]*c.y[i]);
	printf("%12.4f %12.4f %12.4f %12.4f %10.4f %10.4f\n",c.x[i],c.y[i],c.rx[i],c.ry[i],c.xres[i],c.yres[i]);
	c.res[i]=sqrt(c.xres[i]*c.xres[i]+c.yres[i]*c.yres[i]);
	c.xrms+=c.xres[i]*c.xres[i];
	c.yrms+=c.yres[i]*c.yres[i];
	c.rms+=c.xres[i]*c.xres[i]+c.yres[i]*c.yres[i];
	m++;
      }
    }
    c.xrms=sqrt(c.xrms/(float) m);
    c.yrms=sqrt(c.yrms/(float) m);
    c.rms=sqrt(c.rms/(float) m);
    
    // Deselect outliers
    for (i=0;i<c.n;i++) {
      if (c.res[i]>2*c.rms)
	c.usage[i]=0;
    }
  }
  printf("%12.8lf %10.6lf %10.6lf %8.4f %8.4f %8.4f %8.4f\n",img.mjd,t.ra0,t.de0,t.a[1],t.a[2],t.b[1],t.b[2]);
  printf("%d/%d %f %f %f\n",m,c.n,c.xrms,c.yrms,c.rms);

  //  add_fits_keywords(t,"test.fits");
  modify_fits_keywords(t,filename);

  return 0;
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

// Read catalog
struct catalog read_catalog(char *filename)
{
  int i=0;
  char line[LIM];
  FILE *file;
  struct catalog c;

  file=fopen(filename,"r");
  while (fgetline(file,line,LIM)>0) {
    sscanf(line,"%f %f %lf %lf",&c.x[i],&c.y[i],&c.ra[i],&c.de[i]);
    c.usage[i]=1;
    
    i++;
  }
  fclose(file);
  c.n=i;

  return c;
}

// Get a x and y from a RA and Decl
void forward(double ra0,double de0,double ra,double de,float *x,float *y)
{
  int i;
  char pcode[4]="TAN";
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;
  double rx,ry;

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
    if (celfwd(pcode,ra,de,&cel,&phi,&theta,&prj,&rx,&ry)) {
      printf("Error in Projection (celfwd)\n");
      return;
    }
  }
  *x=rx*3600.;
  *y=ry*3600.;

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
void reverse(double ra0,double de0,float x,float y,double *ra,double *de)
{
  int i;
  char pcode[4]="TAN";
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;
  double rx,ry;

  rx=x/3600.;
  ry=y/3600.;

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
    if (celrev(pcode,rx,ry,&prj,&phi,&theta,&cel,ra,de)) {
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

// Read fits image
struct image read_fits(char *filename)
{
  int i,j,k,l,m;
  qfitsloader ql;
  char key[FITS_LINESZ+1];
  char val[FITS_LINESZ+1];
  struct image img;

  // Image size
  img.naxis1=atoi(qfits_query_hdr(filename,"NAXIS1"));
  img.naxis2=atoi(qfits_query_hdr(filename,"NAXIS2"));

  // MJD
  //  img.mjd=(double) atof(qfits_query_hdr(filename,"MJD-OBS"));
  img.mjd=0.0;


  return img;
}
