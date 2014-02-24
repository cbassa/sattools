#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <qfits.h>
#include <jpeglib.h>
#include <libexif/exif-data.h>
#include <getopt.h>

struct image {
  int nx,ny,nz;
  float *z;
  double mjd;
  char nfd[32],observer[32];
  int cospar;
  float exptime;
};
struct image read_jpg(char *filename);
void write_fits(struct image img,char *filename);
double date2mjd(int year,int month,double day);
double nfd2mjd(char *date);
void mjd2nfd(double mjd,char *nfd);

// Read fits image
struct image read_fits(char *filename)
{
  int i,j,k,l,m;
  qfitsloader ql;
  char key[FITS_LINESZ+1] ;
  struct image img;
  float s1,s2,avg,std;

  // Set plane
  ql.xtnum = 0;
  ql.pnum = 0;

  // Set loadtype
  ql.ptype = PTYPE_FLOAT;

  // Set filename
  ql.filename=filename;

  // Image size
  img.nx=atoi(qfits_query_hdr(filename,"NAXIS1"));
  img.ny=atoi(qfits_query_hdr(filename,"NAXIS2"));
  img.nz=1;

  // Initialize load
  if (qfitsloader_init(&ql) != 0) 
    printf("Error initializing data loading\n");

  // Test load
  if (qfits_loadpix(&ql) != 0) 
    printf("Error loading actual data\n");

  // Allocate image memory
  img.z=(float *) malloc(sizeof(float) * img.nx*img.ny*img.nz);

  // Fill z array
  for (i=0,l=0,m=0;i<img.nx;i++) {
    for (j=0;j<img.ny;j++) {
      img.z[l]=ql.fbuf[l];
      l++;
    }
  }

  return img;
}


struct image rebin(struct image raw,int nbin)
{
  int i,j,k;
  int ii,jj,kk;
  struct image img;

  img.nx=raw.nx/nbin;
  img.ny=raw.ny/nbin;
  img.nz=1;
  img.z=(float *) malloc(sizeof(float)*img.nx*img.ny*img.nz);

  for (i=0;i<img.nx;i++) {
    for (j=0;j<img.ny;j++) {
      k=i+img.nx*j;
      img.z[k]=0.0;
      for (ii=0;ii<nbin;ii++) {
	for (jj=0;jj<nbin;jj++) {
	  kk=ii+nbin*i+raw.nx*(jj+nbin*j);
	  img.z[k]+=raw.z[kk];
	}
      }
    }
  }
  img.mjd=raw.mjd;
  img.cospar=raw.cospar;
  img.exptime=raw.exptime;
  strcpy(img.nfd,raw.nfd);
  strcpy(img.observer,raw.observer);
    

  return img;
}

void usage(void)
{
  printf("jpg2fits i:t:o:d:Z:c:T:O:b:h\n\n");
  printf("-i   input JPG file\n");
  printf("-o   output FITS file\n");
  printf("-t   start time (YYYY-MM-DDTHH:MM:SS.sss)\n");
  printf("-d   delay (in seconds)\n");
  printf("-Z   timezone offset (in seconds)\n");
  printf("-T   exposure time (in seconds)\n");
  printf("-c   COSPAR site number\n");
  printf("-O   observer name\n");
  printf("-b   binning factor\n");
  printf("-h   this help\n");

  exit(0);

  return;
}

int main(int argc,char *argv[])
{
  int arg;
  struct image img,raw;
  char infile[64],outfile[64]="",nfd[32]="2000-01-01T00:00:00";
  double mjd=51544.0,delay=0.0,tz=0.0;
  int cospar=0;
  char observer[32]="Cees Bassa";
  float exptime=10.06;
  int flag=0,nbin=1,readfits=0;

  // Decode options
  if (argc>1) {
    while ((arg=getopt(argc,argv,"i:t:o:d:Z:c:T:O:b:h"))!=-1) {
      switch(arg) {
	
      case 'i':
	strcpy(infile,optarg);
	break;
	
      case 'd':
	delay=(double) atof(optarg);
	break;
	
      case 'b':
	nbin=atoi(optarg);
	break;
	
      case 'Z':
	tz=(double) atof(optarg);
	break;
	
      case 'o':
	strcpy(outfile,optarg);
	flag=1;
	break;
	
      case 't':
	strcpy(nfd,optarg);
	break;
	
      case 'c':
	cospar=atoi(optarg);
	break;
	
      case 'T':
	exptime=atof(optarg);
	break;
	
      case 'O':
	strcpy(observer,optarg);
	break;
	
      case 'h':
	usage();
	break;
	
      default:
	usage();
	return 0;
	
      }
    }
  } else {
    usage();
  }

  if (infile!=NULL) {
    if (nbin==1) {
      if (readfits==0)
	img=read_jpg(infile);
      else
	img=read_fits(infile);
    } else {
      if (readfits==0)
	raw=read_jpg(infile);
      else
	raw=read_fits(infile);
      img=rebin(raw,nbin);
    }
  }

  if (nfd!=NULL) {
    // Compute time
    mjd=nfd2mjd(nfd);
    mjd+=(delay+tz)/86400.0;
    mjd2nfd(mjd,nfd);

    // Into file
    strcpy(img.nfd,nfd);
    img.mjd=mjd;
  }

  if (flag==0)
    sprintf(outfile,"%s.fits",img.nfd);

  // Set properties
  img.cospar=cospar;
  img.exptime=exptime;
  strcpy(img.observer,observer);

  if (outfile!=NULL)
    write_fits(img,outfile);

  // Free
  free(img.z);
  if (nbin!=1)
    free(raw.z);

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
  ExifData *ed;
  ExifEntry *entry;

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
      k=i+(img.ny-j-1)*img.nx;
      img.z[k]=0.0;
      for (l=0;l<img.nz;l++) {
	m=img.nz*(i+img.nx*j)+l;
	img.z[k]+=(float) raw_image[m];
      }
      img.z[k]/=3.0;
    }
  }

  // Free allocated memory
  free(row_pointer[0]);
  free(raw_image);

  // Close file
  fclose(file);

  /*
  // Get exif info
  ed=exif_data_new_from_file(filename);
  if (!ed) {
    printf("File not readable or no EXIF data in file %s\n",filename);
  } else {
    entry=exif_content_get_entry(ed->ifd[0],EXIF_TAG_DATE_TIME);
    exif_entry_get_value(entry,img.nfd, sizeof(img.nfd));
    img.nfd[4]='-';
    img.nfd[7]='-';
    img.nfd[10]='T';
    img.nfd[20]='\0';
    img.mjd=nfd2mjd(img.nfd);
  }
  */
  return img;
}

// Write fits file
void write_fits(struct image img,char *filename)
{
  int i,j,k,l;
  int *ibuf;
  qfitsdumper qd;
  qfits_header *qh;
  char key[FITS_LINESZ+1] ;
  char val[FITS_LINESZ+1] ;
  char com[FITS_LINESZ+1] ;
  char lin[FITS_LINESZ+1] ;
  FILE *file;

  // Create FITS header
  qh=qfits_header_default();

  // Add stuff
  qfits_header_add(qh,"BITPIX","16"," ",NULL);
  qfits_header_add(qh,"NAXIS","2"," ",NULL);
  sprintf(val,"%i",img.nx);
  qfits_header_add(qh,"NAXIS1",val," ",NULL);
  sprintf(val,"%i",img.ny);
  qfits_header_add(qh,"NAXIS2",val," ",NULL);
  qfits_header_add(qh,"BSCALE","1.0"," ",NULL);
  qfits_header_add(qh,"BZERO","0.0"," ",NULL);
  qfits_header_add(qh,"DATAMAX","255.0"," ",NULL);
  qfits_header_add(qh,"DATAMIN","0.0"," ",NULL);

  // Astrometry keywors
  sprintf(val,"%f",img.nx/2.0);
  qfits_header_add(qh,"CRPIX1",val," ",NULL);
  sprintf(val,"%f",img.ny/2.0);
  qfits_header_add(qh,"CRPIX2",val," ",NULL);
  qfits_header_add(qh,"CRVAL1","0.0"," ",NULL);
  qfits_header_add(qh,"CRVAL2","0.0"," ",NULL);
  qfits_header_add(qh,"CD1_1","0.0"," ",NULL);
  qfits_header_add(qh,"CD1_2","0.0"," ",NULL);
  qfits_header_add(qh,"CD2_1","0.0"," ",NULL);
  qfits_header_add(qh,"CD2_2","0.0"," ",NULL);
  qfits_header_add(qh,"CTYPE1","'RA---TAN'"," ",NULL);
  qfits_header_add(qh,"CTYPE2","'DEC--TAN'"," ",NULL);
  qfits_header_add(qh,"CUNIT1","'deg'"," ",NULL);
  qfits_header_add(qh,"CUNIT2","'deg'"," ",NULL);
  qfits_header_add(qh,"CRRES1","0.0"," ",NULL);
  qfits_header_add(qh,"CRRES2","0.0"," ",NULL);
  qfits_header_add(qh,"EQUINOX","2000.0"," ",NULL);
  qfits_header_add(qh,"RADECSYS","ICRS"," ",NULL);
  sprintf(val,"%s",img.nfd);
  qfits_header_add(qh,"DATE-OBS",val," ",NULL);
  sprintf(val,"%lf",img.mjd);
  qfits_header_add(qh,"MJD-OBS",val," ",NULL);
  sprintf(val,"%d",img.cospar);
  qfits_header_add(qh,"COSPAR",val," ",NULL);
  sprintf(val,"%f",img.exptime);
  qfits_header_add(qh,"EXPTIME",val," ",NULL);
  sprintf(val,"%s",img.observer);
  qfits_header_add(qh,"OBSERVER",val," ",NULL);

  // Dump fitsheader
  //  qfits_header_dump(qh,stdout);
  
  // Dump to file
  file=fopen(filename,"w");
  qfits_header_dump(qh,file);
  fclose(file);


  // Fill buffer
  ibuf=malloc(img.nx*img.ny*sizeof(int));
  for (i=0,l=0;i<img.nx;i++) {
    for (j=img.ny-1;j>=0;j--) {
      ibuf[l]=(int) img.z[l];
      
      l++;
    }
  }

  // Set parameters
  qd.filename=filename;
  qd.npix=img.nx*img.ny;
  qd.ptype=PTYPE_INT;
  qd.ibuf=ibuf;
  qd.out_ptype=BPP_16_SIGNED;

  // Dump
  qfits_pixdump(&qd);

  free(ibuf);

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

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min;
  float sec;
  double mjd,dday;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%f",&year,&month,&day,&hour,&min,&sec);

  dday=day+hour/24.0+min/1440.0+sec/86400.0;
  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Compute Date from Julian Day
void mjd2nfd(double mjd,char *nfd)
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

  sprintf(nfd,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);

  return;
}
