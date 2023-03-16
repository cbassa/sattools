#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qfits.h"
#include <cpgplot.h>
#include <wcslib/cel.h>
#include <jpeglib.h>
#include <getopt.h>

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
};
struct jpeg_image {
  int nx,ny,nz;
  float *z;
};
struct jpeg_image read_jpg(char *filename);
void write_jpg(char *filename,struct jpeg_image img);
struct image read_fits(char *filename,int pnum);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double ra0,double de0,double x,double y,double *ra,double *de);

int main(int argc,char *argv[])
{
  int i,j,k,l,m;
  struct image img;
  struct jpeg_image jpg,out;
  double rx,ry,ra,de,rx0,ry0;
  double x,y,d;
  double drx=-10.0,dry=10.0;
  double ra0=237.0,de0=12.5;
  int arg=0;
  char *fitsfile,*jpgfile,*outfile;

    // Decode options
  while ((arg=getopt(argc,argv,"j:f:o:R:D:s:"))!=-1) {
    switch(arg) {

    case 'j':
      jpgfile=optarg;
      break;

    case 'f':
      fitsfile=optarg;
      break;

    case 'o':
      outfile=optarg;
      break;

    case 'R':
      ra0=atof(optarg);
      break;

    case 'D':
      de0=atof(optarg);
      break;

    case 's':
      dry=atof(optarg);
      drx=-dry;
      break;

    default:
      return 0;
    }
  }


  // Read image
  img=read_fits(fitsfile,0);
  jpg=read_jpg(jpgfile);

  out.nx=3000;
  out.ny=6000;
  out.nz=3;

  /*
  img.x0*=4.0;
  img.y0*=4.0;
  img.a[1]/=4.0;
  img.a[2]/=4.0;
  img.b[1]/=4.0;
  img.b[2]/=4.0;
  */
  out.z=(float *) malloc(sizeof(float)*out.nx*out.ny*out.nz);

  for (i=0;i<out.nx;i++) {
    for (j=0;j<out.ny;j++) {
      // Set rx,ry
      rx=drx*(float) (i-0.5*out.nx);
      ry=dry*(float) (j-0.5*out.ny);

      // Obtain ra/dec for output image
      reverse(ra0,de0,rx,ry,&ra,&de);
      
      // Obtain rx/ry for input image
      forward(img.ra0,img.de0,ra,de,&rx0,&ry0);

      // Compute pixel position
      d=img.a[1]*img.b[2]-img.a[2]*img.b[1];
      x=(+rx0*img.b[2]-ry0*img.a[2])/d+img.x0;
      y=(-rx0*img.b[1]+ry0*img.a[1])/d+img.y0;

      // Fill image
      for (k=0;k<jpg.nz;k++) {
	l=out.nz*(i+out.nx*(out.ny-j-1))+k;
	m=jpg.nz*((int) x+jpg.nx*(int) (jpg.ny-y-1))+k;
	if (x>0.0 && x<jpg.nx && y>0.0 && y<jpg.ny) 
	  out.z[l]=jpg.z[m];
	else
	  out.z[l]=0.0;
	
      }
    }
  }

  // Dump
  write_jpg(outfile,out);

  // Free
  free(img.z);
  free(jpg.z);
  free(out.z);

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

struct jpeg_image read_jpg(char *filename)
{
  int i=0,j,k,l,m;
  unsigned long location=0;
  struct jpeg_image img;
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];
  unsigned char *raw_image=NULL;
  FILE *file;

  // Open file
  file=fopen(filename,"rb");
  if (!file)
    perror("Error opening file");

  // Get header info, decompress
  cinfo.err=jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo,file);
  jpeg_read_header(&cinfo,TRUE);
  jpeg_start_decompress(&cinfo);

  // Allocate memory
  raw_image=(unsigned char *) malloc(cinfo.output_width*cinfo.output_height*cinfo.num_components);

  // Read image, one scan at a time
  row_pointer[0]=(unsigned char *) malloc(cinfo.output_width*cinfo.num_components);
  while(cinfo.output_scanline<cinfo.image_height) {
    jpeg_read_scanlines(&cinfo,row_pointer,1);
    for(i=0;i<cinfo.image_width*cinfo.num_components;i++) 
      raw_image[location++]=row_pointer[0][i];
  }
  // wrap up decompression, destroy objects, free pointers and close open files
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Copy image to image struct
  img.nx=cinfo.image_width;
  img.ny=cinfo.image_height;
  img.nz=cinfo.num_components;
  img.z=(float *) malloc(sizeof(float)*img.nx*img.ny*img.nz);

  // Fill image
  for (i=0;i<img.nx;i++) {
    for (j=0;j<img.ny;j++) {
      for (k=0;k<img.nz;k++) {
	l=img.nz*(i+img.nx*j)+k;
	img.z[l]=(float) raw_image[l];
      }
    }
  }

  // Free allocated memory
  free(row_pointer[0]);
  free(raw_image);

  // Close file
  fclose(file);

  return img;
}

// Write jpg
void write_jpg(char *filename,struct jpeg_image img)
{
  int i,j,k,l,m;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];
  FILE *outfile;
  unsigned char *raw_image=NULL;

  outfile=fopen(filename,"wb");
  cinfo.err=jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo,outfile);
  cinfo.image_width=img.nx;	
  cinfo.image_height=img.ny;
  cinfo.input_components=3;
  cinfo.in_color_space=JCS_RGB;
  jpeg_set_defaults(&cinfo);
  jpeg_start_compress(&cinfo,TRUE);

  // Allocate memory
  raw_image=(unsigned char *) malloc(cinfo.image_width*cinfo.image_height*cinfo.input_components);

  // Fill image
  for (i=0;i<img.nx;i++) {
    for (j=0;j<img.ny;j++) {
      for (k=0;k<img.nz;k++) {
	l=img.nz*(i+img.nx*j)+k;
	raw_image[l]=(unsigned char) img.z[l];
      }
    }
  }

  while(cinfo.next_scanline<cinfo.image_height) {
    row_pointer[0]=&raw_image[cinfo.next_scanline*cinfo.image_width*cinfo.input_components];
    jpeg_write_scanlines(&cinfo,row_pointer,1);
  }
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
  fclose(outfile);

  return;
}

