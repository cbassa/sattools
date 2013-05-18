#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cel.h"
#include "cpgplot.h"
#include "qfits.h"

#define LIM 80
#define NMAX 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI

struct image {
  char filename[64];
  int naxis1,naxis2,naxis3,nframes;
  float *zavg,*zstd,*zmax,*znum,*ztrk;
  double ra0,de0;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  double mjd;
  float *dt,exptime;
  char nfd[32];
  int cospar;
};
struct image read_fits(char *filename);
void write_pgm(char *filename,struct image img);

int main(int argc,char *argv[])
{
  int i;
  struct image img;

  img=read_fits(argv[1]);

  write_pgm("avg.pgm",img);
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
  img.naxis3=atoi(qfits_query_hdr(filename,"NAXIS3"));

  // MJD
  img.mjd=(double) atof(qfits_query_hdr(filename,"MJD-OBS"));
  strcpy(img.nfd,qfits_query_hdr(filename,"DATE-OBS"));

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
  img.exptime=atof(qfits_query_hdr(filename,"EXPTIME"));
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
  if (img.naxis3==5) 
    img.ztrk=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);

  // Set parameters
  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;

  // Loop over planes
  for (k=0;k<img.naxis3;k++) {
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
	if (k==1) img.zstd[l]=ql.fbuf[l];
	if (k==2) img.zmax[l]=ql.fbuf[l];
	if (k==3) img.znum[l]=ql.fbuf[l];
	if (img.naxis3==5) {
	  if (k==0) img.ztrk[l]=ql.fbuf[l];
	  if (k==4) img.zavg[l]=ql.fbuf[l];
	} else {
	  if (k==0) img.zavg[l]=ql.fbuf[l];
	}

	l++;
      }
    }
  }

  return img;
}

// Write pgm file
void write_pgm(char *filename,struct image img)
{
  int i,j,k,l,n;
  FILE *file;
  float z,zavgmin,zavgmax,zstdmin,zstdmax,zmaxmin,zmaxmax;
  float s1,s2,avg,std;
  unsigned char *buffer;

  n=img.naxis1*img.naxis2;
  for (j=0;j<3;j++) {
    for (i=0,s1=0.0,s2=0.0;i<n;i++) {
      if (j==0) z=img.zavg[i];
      if (j==1) z=img.zstd[i];
      if (j==2) z=img.zmax[i];
      s1+=z;
      s2+=z*z;
    }
    avg=s1/(float) n;
    std=sqrt(s2/(float) n-avg*avg);
    if (j==0) {
      zavgmin=avg-2*std;
      zavgmax=avg+3*std;
    }
    if (j==1) {
      zstdmin=avg-2*std;
      zstdmax=avg+3*std;
    }
    if (j==2) {
      zmaxmin=avg-2*std;
      zmaxmax=avg+3*std;
    }
  }
  buffer=(unsigned char *) malloc(sizeof(unsigned char)*4*n);

  for (j=0,l=0;j<img.naxis2;j++) {
    for (i=0;i<img.naxis1;i++) {
      k=i+(img.naxis2-j-1)*img.naxis1;
      z=255.0*(img.zavg[k]-zavgmin)/(zavgmax-zavgmin);
      if (z>=255.0)
	z=255.0;
      if (z<0.0)
	z=0.0;
      buffer[l++]=(unsigned char) z;
    }
    for (i=0;i<img.naxis1;i++) {
      k=i+(img.naxis2-j-1)*img.naxis1;
      z=255.0*(img.zstd[k]-zstdmin)/(zstdmax-zstdmin);
      if (z>=255.0)
	z=255.0;
      if (z<0.0)
	z=0.0;
      buffer[l++]=(unsigned char) z;
    }
  }
  for (j=0;j<img.naxis2;j++) {
    for (i=0;i<img.naxis1;i++) {
      k=i+(img.naxis2-j-1)*img.naxis1;
      z=255*(img.zmax[k]-zmaxmin)/(zmaxmax-zmaxmin);
      if (z>=255.0)
	z=255.0;
      if (z<0.0)
	z=0.0;
      buffer[l++]=(unsigned char) z;
    }
    for (i=0;i<img.naxis1;i++) {
      k=i+(img.naxis2-j-1)*img.naxis1;
      z=img.znum[k];
      if (z>=255.0)
	z=255.0;
      if (z<0.0)
	z=0.0;
      buffer[l++]=(unsigned char) z;
    }
  }
  file=fopen(filename,"wb");
  fprintf(file,"P5\n%d %d\n255\n",2*img.naxis1,2*img.naxis2);
  fwrite(buffer,4*n,sizeof(unsigned char),file);
  fclose(file);

  return;
}

// Write pgm file
void write_pgm2(char *filename,struct image img)
{
  int i,j,k;
  FILE *file;
  float z;

  file=fopen(filename,"w");
  fprintf(file,"P5\n# %.23s\n%d %d\n255\n",img.nfd+1,img.naxis1,img.naxis2);
  for (j=0;j<img.naxis2;j++) {
    for (i=0;i<img.naxis1;i++) {
      k=i+(img.naxis2-j-1)*img.naxis1;
      //      z=255.0*(img.zavg[k]-30.0)/(60.0-30.0);
      z=img.zstd[k];
      //z=255.0*(img.ztrk[k]-30.0)/(60.0-30.0);
      if (z>255.0)
	z=255.0;
      if (z<0.0)
	z=0.0;
      fprintf(file,"%c",(char) z);
    }
  }
  fclose(file);

  return;
}
