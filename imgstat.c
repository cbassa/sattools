#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qfits.h"

struct image {
  char filename[64];
  int naxis1,naxis2,nframes;
  float *zavg,*zstd,*zmax,*znum;
  double ra0,de0;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  double mjd;
  float *dt,exptime;
  char nfd[32];
  int cospar;
};
struct image read_fits(char *filename);

int main(int argc,char *argv[])
{
  int i;
  struct image img;
  float zavg,zstd,sx,sy,wx,wy;

  // Read image
  img=read_fits(argv[1]);

  // Compute statistics
  for (i=0,zavg=0.0;i<img.naxis1*img.naxis2;i++)
    zavg+=img.zmax[i];
  zavg/=(float) img.naxis1*img.naxis2;
  for (i=0,zstd=0.0;i<img.naxis1*img.naxis2;i++)
    zstd+=pow(img.zmax[i]-zavg,2);
  zstd=sqrt(zstd/(float) (img.naxis1*img.naxis2));

  // Image scale
  sx=sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1]);
  sy=sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2]);
  wx=img.naxis1*sx/3600.0;
  wy=img.naxis2*sy/3600.0;

  printf("%s %14.8lf %10.6f %10.6f %.2f %.2f %.2f %.2f %.2f %.2f %.1f %.1f\n",argv[1],img.mjd,img.ra0,img.de0,wx,wy,sx,sy,img.xrms,img.yrms,zavg,zstd);

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

  // Set parameters
  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;

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

  return img;
}

