#include <stdio.h>
#include <stdlib.h>
#include <jpeglib.h>
#include <math.h>
#include <cpgplot.h>

#define NMAX 1024

struct image {
  int nx,ny,nz;
  float *z;
};
struct image read_jpg(char *filename);
void write_jpg(char *filename,struct image img);


int main(int argc,char *argv[])
{
  int i,j,flag,n;
  struct image avg,max,raw;
  
  for (flag=0,i=1,n=0;i<argc;i++,n++) {
    printf("%d %s\n",i,argv[i]);
    // Read image
    raw=read_jpg(argv[i]);    

    // If first image, initialize
    if (flag==0) {
      avg.nx=raw.nx;
      avg.ny=raw.ny;
      avg.nz=raw.nz;
      avg.z=(float *) malloc(sizeof(float)*avg.nx*avg.ny*avg.nz);
      for (j=0;j<avg.nx*avg.ny*avg.nz;j++)
	avg.z[j]=0.0;

      max.nx=raw.nx;
      max.ny=raw.ny;
      max.nz=raw.nz;
      max.z=(float *) malloc(sizeof(float)*max.nx*max.ny*max.nz);
      for (j=0;j<max.nx*max.ny*max.nz;j++)
	max.z[j]=0.0;

      flag=1;
    }
    
    // Add values
    for (j=0;j<avg.nx*avg.ny*avg.nz;j++)
      avg.z[j]+=raw.z[j];

    // Add values
    for (j=0;j<avg.nx*avg.ny*avg.nz;j++)
      if (raw.z[j]>max.z[j])
	max.z[j]=raw.z[j];

    // Free
    free(raw.z);
  }

  // Average;
  for (j=0;j<avg.nx*avg.ny*avg.nz;j++)
    avg.z[j]/=(float) n;

  // Write
  write_jpg("average.jpg",avg);
  write_jpg("maximum.jpg",max);

  // Free
  free(avg.z);
  free(max.z);

  return 0;
}

struct image read_jpg(char *filename)
{
  int i=0,j,k,l,m;
  unsigned long location=0;
  struct image img;
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
void write_jpg(char *filename,struct image img)
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
