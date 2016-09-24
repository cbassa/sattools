#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <qfits.h>
#include <getopt.h>

struct image {
  int nx,ny;
  char timestamp[64];
  unsigned char *c;
};
struct fourframe {
  int nx,ny,nt,nlayer;
  char timestamp[64],observer[64];
  double mjd;
  float *z,*dt;
  int cospar;
};
int fgetline(FILE *file,char *s,int lim);
void write_fits(char *filename,struct fourframe ff);
struct image read_pgm(char *filename,int *status);
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void write_pgm(char *filename,struct fourframe ff);
void mjd2date(double mjd,char *date);

void usage(void)
{
  printf("pgm2fits p:w:h:s:n:Dd:x:y:c:o:gm:t:r:\n\n");
  printf("-p   image prefix\n");
  printf("-w   image width in pixels\n");
  printf("-h   image height in pixels\n");
  printf("-s   number of first image to process\n");
  printf("-n   number of images to process\n");
  printf("-D   toggle for creating dark frame\n");
  printf("-d   toggle for creating dark frame\n");
  printf("-d   filename of dark frame to substract\n");
  printf("-m   filename of mask frame to apply\n");
  printf("-x   tracking rate in x (pix/s)\n");
  printf("-y   tracking rate in y (pix/s)\n");
  printf("-c   COSPAR [default 4553]\n");
  printf("-o   observer [default \"Cees Bassa\"]\n");
  printf("-g   toggle for guiding??\n");
  printf("-t   time stamp of first image [YYYY-MM-DDTHH:MM:SS.SSS]\n");
  printf("-r   frame rate (frames/s)\n");
  exit(0);
}

int main(int argc,char *argv[])
{
  int i,j,k,l,m,k0=1,status,nt,darkout=0,darkin=0,track=0,maskin=0;
  int n,n0,di,dj,npix;
  struct image *img,drk,msk;
  struct fourframe ff;
  char filename[128],nfd[64];
  float s1,s2,z;
  float avg,std,max,cnt,*trk;
  int *wt;
  int arg=0;
  char *path,*darkfile,*maskfile;
  double mjd,mjd0=0.0;
  float dxdn=0.0,dydn=0.0,dx,dy;
  int guide=0,timereset=0;
  float framerate=25.0;
  char *env;

  // Set defaults
  env=getenv("ST_COSPAR");
  ff.cospar=atoi(env);
  strcpy(ff.observer,"Cees Bassa");

  // Decode options
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:w:h:s:n:Dd:x:y:c:o:gm:t:r:"))!=-1) {
      switch(arg) {
      case 'p':
	path=optarg;
	break;

      case 'g':
	guide=1;
	break;

      case 'w':
	ff.nx=atoi(optarg);
	break;
	
      case 'h':
	ff.ny=atoi(optarg);
	break;

      case 'c':
	ff.cospar=atoi(optarg);
	break;

      case 'o':
	strcpy(ff.observer,optarg);
	break;

      case 's':
	k0=atoi(optarg);
	break;
	
      case 'n':
	nt=atoi(optarg);
	break;

      case 'D':
	darkout=1;
	break;
	
      case 'd':
	darkin=1;
	darkfile=optarg;
	break;

      case 'm':
	maskin=1;
	maskfile=optarg;
	break;
	
      case 'x':
	dxdn=atof(optarg);
	track=1;
	break;
	
      case 'y':
	dydn=atof(optarg);
	track=1;
	break;

      case 't':
	strcpy(nfd,optarg);
	mjd0=nfd2mjd(nfd);
	timereset=1;
	break;

      case 'r':
	framerate=atof(optarg);
	timereset=1;
	break;

      default:
	usage();
      }
    }
  } else {
    usage();
  }

  // Add layer
  if (track==1) 
    ff.nlayer=5;
  else
    ff.nlayer=4;

  // Allocate
  ff.z=(float *) malloc(ff.nlayer*sizeof(float)*ff.nx*ff.ny);
  ff.dt=(float *) malloc(sizeof(float)*nt);
  trk=(float *) malloc(sizeof(float)*ff.nx*ff.ny);
  wt=(int *) malloc(sizeof(float)*ff.nx*ff.ny);
  img=(struct image *) malloc(sizeof(struct image)*nt);

  // Read dark file
  if (darkin==1) 
    drk=read_pgm(darkfile,&status);

  // Read mask file
  if (maskin==1)
    msk=read_pgm(maskfile,&status);

  // Loop over files
  for (k=0,l=0;k<nt;k++) {
    sprintf(filename,"%s%06d.pgm",path,k+k0);
    img[l]=read_pgm(filename,&status);

    // Reset time
    if (timereset==1) {
      mjd=mjd0+(double) (k+k0)/(86400.0*framerate);
      mjd2date(mjd,img[l].timestamp);
    }

    ff.dt[l]=86400.0*(nfd2mjd(img[l].timestamp)-nfd2mjd(img[0].timestamp));
    if (status==0) {
      printf("Read %s\n",filename);
      l++;
    } else {
      break;
    }
  }
  ff.nt=l;
  strcpy(ff.timestamp,img[0].timestamp);
  ff.mjd=nfd2mjd(img[0].timestamp);

  printf("Accumulating image statistics\n");
  // Loop over pixels
  for (i=0;i<ff.nx;i++) {
    for (j=0;j<ff.ny;j++) {
      n=i+ff.nx*j;
      
      s1=0.0;
      s2=0.0;
      max=0.0;
      cnt=0.0;

      // Loop over images
      for (k=0;k<ff.nt;k++) {
	if (darkin==0)
	  z=(float) img[k].c[n];
	else if (darkin==1)
	  z=(float) (img[k].c[n]-drk.c[n]);
	s1+=z;
	s2+=z*z;
	if (z>max) {
	  max=z;
	  cnt=(float) k;
	}
      }
      s1-=max;
      s2-=max*max;
      avg=s1/(float) ff.nt;
      std=sqrt((s2-s1*avg)/(float) (ff.nt-1));
      
      // Reset masked pixels
      if (maskin==1 && msk.c[n]==0.0) {
	avg=0.0;
	std=0.0;
	max=0.0;
	cnt=128.0;
      }

      for (m=0;m<ff.nlayer;m++) {
	l=i+(ff.ny-j-1)*ff.nx+m*ff.nx*ff.ny;
	if (m==0) ff.z[l]=avg;
	if (m==1) ff.z[l]=std;
	if (m==2) ff.z[l]=max;
	if (m==3) ff.z[l]=cnt;
	if (m==4) ff.z[l]=avg;
      }
    }
  }

  // Create tracked layer
  if (track==1) {
    printf("Creating tracked layer\n");

    // Set weights
    for (i=0;i<ff.nx*ff.ny;i++) {
      wt[i]=0;
      trk[i]=0.0;
    }

    // Loop over frames
    for (l=0;l<ff.nt;l++) {
      // Offset
      dx=dxdn*(l-ff.nt/2);
      dy=dydn*(l-ff.nt/2);
      
      // Integer offset
      di=(int) floor(dx+0.5);
      dj=(int) floor(dy+0.5);

      // Set
      for (i=0;i<ff.nx;i++) {
	for (j=0;j<ff.ny;j++) {
	  k=i+ff.nx*j;
	  k0=i+di+ff.nx*(j+dj);
	  if (i+di>0 && i+di<ff.nx && j+dj>0 && j+dj<ff.ny) {
	    wt[k]+=1;
	    trk[k]+=(float) img[l].c[k0];
	  }
	}
      }
    }
    // Save layer
    for (i=0;i<ff.nx;i++) {
      for (j=0;j<ff.ny;j++) {
	k=i+ff.nx*j;
	if (guide==0)
	  l=i+(ff.ny-j-1)*ff.nx;
	else
	  l=i+(ff.ny-j-1)*ff.nx+4*ff.nx*ff.ny;
	if (wt[k]>0)
	  ff.z[l]=trk[k]/(float) wt[k];
	else
	  ff.z[l]=trk[k];
      }
    }
  }

  // Write fits
  if (ff.timestamp!=NULL)
    sprintf(filename,"%s.fits",ff.timestamp);
  else
    strcpy(filename,"test.fits");
  write_fits(filename,ff);

  // Write dark frame
  if (darkout==1) 
    write_pgm("dark.pgm",ff);

  // Free
  for (k=0;k<ff.nt;k++)
    free(img[k].c);
  free(trk);
  free(img);
  free(wt);
  free(ff.dt);
  free(ff.z);

  return 0;
}

// Read pgm file
struct image read_pgm(char *filename,int *status)
{
  int i;
  struct image img;
  FILE *file;
  char hbuf[64];

  // Open file
  file=fopen(filename,"rb");
  if (file==NULL) {
    *status=1;
    return img;
  }

  // Read PGM format
  fgetline(file,hbuf,64);
  if (strcmp(hbuf,"P5")!=0) {
    printf("Not a valid PGM file!\n");
    exit(0);
  }

  // Read timestamp/image size
  fgetline(file,hbuf,64);
  if (strstr(hbuf,"#")!=NULL) {
    strcpy(img.timestamp,hbuf+2);
    fgetline(file,hbuf,64);
    sscanf(hbuf,"%d %d",&img.nx,&img.ny);
  } else {
    strcpy(img.timestamp,"2012-01-01T00:00:00");
    sscanf(hbuf,"%d %d",&img.nx,&img.ny);
  }
  fgetline(file,hbuf,64);

  // Allocate
  img.c=(unsigned char *) malloc(sizeof(unsigned char)*img.nx*img.ny);
  
  // Read
  fread(img.c,1,img.nx*img.ny,file);

  // Close file
  fclose(file);

  *status=0;
  return img;
}

// Get line
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim>0 && (c=fgetc(file))!=EOF && c!='\n')
    s[i++]=c;
  //  if (c=='\n')
  //    s[i++]=c;
  s[i]='\0';

  return i;
}
  
// Write fits file
void write_fits(char *filename,struct fourframe ff)
{
  int i,j,k,l;
  float *fbuf;
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
  qfits_header_add(qh,"BITPIX","-32"," ",NULL);
  //  qfits_header_add(qh,"BITPIX","16"," ",NULL);
  qfits_header_add(qh,"NAXIS","3"," ",NULL);
  sprintf(val,"%i",ff.nx);
  qfits_header_add(qh,"NAXIS1",val," ",NULL);
  sprintf(val,"%i",ff.ny);
  qfits_header_add(qh,"NAXIS2",val," ",NULL);
  sprintf(val,"%i",ff.nlayer);
  qfits_header_add(qh,"NAXIS3",val," ",NULL);
  qfits_header_add(qh,"BSCALE","1.0"," ",NULL);
  qfits_header_add(qh,"BZERO","0.0"," ",NULL);
  qfits_header_add(qh,"DATAMAX","255.0"," ",NULL);
  qfits_header_add(qh,"DATAMIN","0.0"," ",NULL);
  sprintf(val,"%s",ff.timestamp);
  qfits_header_add(qh,"DATE-OBS",val," ",NULL);
  // MJD-OBS
  sprintf(val,"%lf",ff.mjd);
  qfits_header_add(qh,"MJD-OBS",val," ",NULL);
  sprintf(val,"%f",ff.dt[ff.nt-1],ff.dt[0]);
  qfits_header_add(qh,"EXPTIME",val," ",NULL);
  sprintf(val,"%d",ff.nt);
  qfits_header_add(qh,"NFRAMES",val," ",NULL);

  // Astrometry keywors
  sprintf(val,"%f",ff.nx/2.0);
  qfits_header_add(qh,"CRPIX1",val," ",NULL);
  sprintf(val,"%f",ff.ny/2.0);
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
  sprintf(val,"%d",ff.cospar);
  qfits_header_add(qh,"COSPAR",val," ",NULL);
  sprintf(val,"'%s'",ff.observer);
  qfits_header_add(qh,"OBSERVER",val," ",NULL);

  // Add timestamps
  for (k=0;k<ff.nt;k++) {
    sprintf(key,"DT%04d",k);
    sprintf(val,"%f",ff.dt[k]);
    qfits_header_add(qh,key,val," ",NULL);
  }

  // Dummy keywords
  for (k=0;k<10;k++) {
    sprintf(key,"DUMY%03d",k);
    qfits_header_add(qh,key,"0.0"," ",NULL);
  }

  // Dump fitsheader
  //  qfits_header_dump(qh,stdout);
  
  // Dump to file
  file=fopen(filename,"w");
  qfits_header_dump(qh,file);
  fclose(file);


  // Fill buffer
  fbuf=malloc(ff.nlayer*ff.nx*ff.ny*sizeof(float));
  //  ibuf=malloc(4*ff.nx*ff.ny*sizeof(int));
  for (i=0,l=0;i<ff.nx;i++) {
    for (j=ff.ny-1;j>=0;j--) {
      for (k=0;k<ff.nlayer;k++) {
	fbuf[l]=ff.z[l];
	//	ibuf[l]=(int) ff.z[l];
   
	l++;
      }
    }
  }

  // Set parameters
  qd.filename=filename;
  qd.npix=ff.nlayer*ff.nx*ff.ny;
  qd.ptype=PTYPE_FLOAT;
  //  qd.ptype=PTYPE_INT;
  qd.fbuf=fbuf;
  //  qd.ibuf=ibuf;
  qd.out_ptype=-32;
  //  qd.out_ptype=BPP_16_SIGNED;

  // Dump
  qfits_pixdump(&qd);

  free(fbuf);
  //  free(ibuf);

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

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min;
  double mjd,dday;
  float sec;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%f",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;
  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Write pgm file
void write_pgm(char *filename,struct fourframe ff)
{
  int i,j,k;
  FILE *file;

  file=fopen(filename,"w");
  fprintf(file,"P5\n# 2013-01-01T00:00:00\n%d %d\n255\n",ff.nx,ff.ny);
  for (j=0;j<ff.ny;j++) {
    for (i=0;i<ff.nx;i++) {
      k=i+(ff.ny-j-1)*ff.nx;
      fprintf(file,"%c",(char) ff.z[k]);
    }
  }
  fclose(file);

  return;
}

// Compute Date from Julian Day
void mjd2date(double mjd,char *date)
{
  double f,jd,dday;
  int z,alpha,a,b,c,d,e;
  int year,month,day,hour,min;
  double sec,x;

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

  sprintf(date,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);
  
  return;
}
