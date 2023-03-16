#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qfits.h"
#include <cpgplot.h>
#include <wcslib/cel.h>
#include <gsl/gsl_multifit.h>
#include <getopt.h>

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
struct image read_fits(char *filename,int pnum);
struct catalog read_pixel_catalog(char *filename,float,float,float);
int fgetline(FILE *file,char *s,int lim);
struct catalog read_astrometric_catalog(char *filename,float mmin,float sx,float sy,float angle);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
int select_nearest(struct catalog c,float x,float y);
void lfit2d(float *x,float *y,float *z,int n,float *a);
void fit_transformation(struct catalog cat,struct catalog ast,int nselect);
struct catalog reread_astrometric_catalog(char *filename,float mmin);
int match_catalogs(struct catalog *cat,struct catalog *ast,float rmax);
double sex2dec(char *s);

void usage(void)
{

  printf("-f <file>      FITS file to calibrate\n");
  printf("-R <hh:mm:ss>  RA of center\n");
  printf("-D <ddd:mm:ss> Decl of center\n");
  printf("-s <scale>     Pixel scale (arcseconds) [10.0]\n");
  printf("-q <parang>    Parallactic angle (degrees) [0.0]\n");
  printf("-m <limmag>    Magnitude limit [9.0]\n");
  printf("-r <radius>    Matching radius (pixels) [10]\n");

  return;
}

int main(int argc,char *argv[])
{
  int i,j;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  int redraw=1,plotcat=0,click=0,nselect=0;
  float x,y,width;
  char c,pixcat[LIM];
  struct catalog cat,ast;
  float sx=-10.0,sy=10.0,q=0.0;
  char *env,starfile[128];
  float r,srmin=1.0,srmax=10.0,smmin=2.0,smmax=8.0,mag=9.0,rmax=10.0,rmask=-1;
  int arg=0;
  char *fitsfile=NULL;
  char sra[16],sde[16];
  double ra,de;
  FILE *file;

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(starfile,"%s/data/tycho2.dat",env);

  // Decode options
  if (argc>1) {
    while ((arg=getopt(argc,argv,"f:R:D:s:q:m:r:M:"))!=-1) {
      switch(arg) {
	
      case 'f':
	fitsfile=optarg;
	break;
	
      case 'R':
	strcpy(sra,optarg);
	if (strchr(sra,':')!=NULL)
	  ra=15.0*sex2dec(sra);
	else
	  ra=atof(sra);
	break;
	
      case 'D':
	strcpy(sde,optarg);
	if (strchr(sde,':')!=NULL)
	  de=sex2dec(sde);
	else
	  de=atof(sde);
	break;
	
      case 's':
	if (strchr(optarg,',')!=NULL) {
	  sscanf(optarg,"%f,%f",&sx,&sy);
	} else {
	  sx=-atof(optarg);
	  sy=-sx;
	}
	break;
	
      case 'q':
	q=atof(optarg);
	break;
	
      case 'm':
	mag=atof(optarg);
	break;
	
      case 'M':
	rmask=atof(optarg);
	break;
	
      case 'r':
	rmax=atof(optarg);
	break;
	
      case 'h':
	usage();
	return 0;
	break;
	
      default:
	usage();
	return 0;
      }
    } 
  } else {
    usage();
    return 0;
  }

  // Read image
  img=read_fits(fitsfile,0);

  // Set variables
  img.ra0=ra;
  img.de0=de;
  img.x0=0.5*(float) img.naxis1;
  img.y0=0.5*(float) img.naxis2;


  // Read pixel catalog
  sprintf(pixcat,"%s.cat",fitsfile);
  cat=read_pixel_catalog(pixcat,img.x0,img.y0,rmask);

  // Read astrometric catalog
  ast=read_astrometric_catalog(starfile,mag,sx,sy,q);

  // Open server
  cpgopen("/xs");
  //  cpgpap(0.,1.0);
  cpgask(0);
  cpgsch(0.8);

  // Default limits
  width=(img.naxis1>img.naxis2) ? img.naxis1 : img.naxis2;
  xmin=0.5*(img.naxis1-width);
  xmax=0.5*(img.naxis1+width);
  ymin=0.5*(img.naxis2-width);
  ymax=0.5*(img.naxis2+width);
  zmin=img.zmin;
  zmax=img.zmax;

  // Forever loop
  for (;;) {
    if (redraw==1) {
      cpgeras();
      
      cpgsvp(0.1,0.95,0.1,0.95);
      cpgwnad(xmin,xmax,ymin,ymax);
      cpglab("x (pix)","y (pix)"," ");
      cpgsfs(2);
      cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);

      cpgimag(img.z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

      // Plot catalog
      if (plotcat==1) {
	cpgsci(3);
	for (i=0;i<cat.n;i++) {
	  if (cat.select[i]!=0) 
	    cpgpt1(cat.x[i],cat.y[i],6);
	  else
	    cpgpt1(cat.x[i],cat.y[i],24);
	}
	cpgsci(1);
      }

      cpgsci(4);
      for (i=0;i<ast.n;i++) {
	r=srmax-(srmax-srmin)*(ast.mag[i]-smmin)/(smmax-smmin);

	// Upscale for image size
	r*=img.naxis1/752.0;
	if (ast.select[i]!=0)
	  cpgpt1(ast.x[i],ast.y[i],6);
	cpgcirc(ast.x[i],ast.y[i],r);
      }
      cpgsci(1);
      redraw=0;
    }

    // Get cursor
    cpgcurs(&x,&y,&c);

    // Quit
    if (c=='q')
      break;

    // Help
    if (c=='h') {
      printf("Calibrates astrometry. Initially requires manual matching of at least three stars. Use 'a' to select star on the image, then 'b' to select star from the catalog (blue circles). If at least three stars are matched, use 'f' and 'm' to converge the calibration.\n\n");
      printf("q     Quit\n");
      printf("a     Select star on image\n");
      printf("b     Select star from catalog\n");
      printf("c     Center image on pixel\n");
      printf("f     Fit calibration\n");
      printf("m     Match stars using current calibration\n");
      printf("z/+   Zoom in on cursor\n");
      printf("x/-   Zoom out on cursor\n");
      printf("p     Plot sextractor catalog\n");
      printf("r     Reset zoom\n");
      
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

    // Center
    if (c=='c') {
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // Fit
    if (c=='f' && nselect>=3) {
      fit_transformation(cat,ast,nselect);
      ast=reread_astrometric_catalog(starfile,mag+1);
      redraw=1;
    }

    // Zoom
    if (c=='z' || c=='+' || c=='=') {
      width/=1.25;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // Match catalogs
    if (c=='m') {
      nselect=match_catalogs(&cat,&ast,rmax);
      redraw=1;
    }

    // Unzoom
    if (c=='x' || c=='-') {
      width*=1.25;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // Plot catalog
    if (c=='p') {
      plotcat = (plotcat==1) ? 0 : 1;
      redraw=1;
    }

    // Reset
    if (c=='r') {
      width=(img.naxis1>img.naxis2) ? img.naxis1 : img.naxis2;
      xmin=0.5*(img.naxis1-width);
      xmax=0.5*(img.naxis1+width);
      ymin=0.5*(img.naxis2-width);
      ymax=0.5*(img.naxis2+width);
      redraw=1;
      continue;
    }

    // Dump
    if (c=='d') {
      file=fopen("dump.dat","w");
      for (i=0;i<nselect;i++) {
	for (j=0;j<cat.n;j++) 
	  if (cat.select[j]==i+1) 
	    fprintf(file,"%lf %lf ",cat.x[j]-img.x0,cat.y[j]-img.y0);
	for (j=0;j<ast.n;j++) 
	  if (ast.select[j]==i+1) 
	    fprintf(file,"%lf %lf\n",ast.ra[j],ast.de[j]);
      }
      fclose(file);
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
  //  img.mjd=atof(qfits_query_hdr(filename,"MJD-OBS"));

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
  for (i=0,s1=0.0;i<img.naxis1*img.naxis2;i++) 
    s1+=img.z[i];
  avg=s1/(float) (img.naxis1*img.naxis2);  
  for (i=0,s2=0.0;i<img.naxis1*img.naxis2;i++) 
    s2+=pow(img.z[i]-avg,2);
  std=sqrt(s2/(float) (img.naxis1*img.naxis2));
  img.zmin=avg-4.0*std;
  img.zmax=avg+16.0*std;
  
  return img;
}

// Read pixel catalog
struct catalog read_pixel_catalog(char *filename,float x0,float y0,float rmin)
{
  int i=0;
  FILE *file;
  char line[LIM];
  struct catalog c;
  float dx,dy,r;

  // Read catalog
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"%s not found!\n",filename);
    exit(0);
  }
  while (fgetline(file,line,LIM)>0) {
    if (i>=NMAX)
      break;
    if (strstr(line,"#")!=NULL) 
      continue;
    sscanf(line,"%f %f %f",&c.x[i],&c.y[i],&c.mag[i]);
    dx=c.x[i]-x0;
    dy=c.y[i]-y0;
    r=sqrt(dx*dx+dy*dy);
    if (r>rmin && rmin>0.0)
      continue;
    c.select[i]=0;
    i++;
  }
  fclose(file);
  c.n=i;

  return c;
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
    forward(img.ra0,img.de0,s.ra,s.de,&rx,&ry);
    x=img.x0+1.0/sx*(cos(angle*D2R)*rx+sin(angle*D2R)*ry);
    y=img.y0+1.0/sy*(-sin(angle*D2R)*rx+cos(angle*D2R)*ry);
    if (x>0.0 && x<(float) img.naxis1 && y>0.0 && y<(float) img.naxis2) {
      c.x[i]=x;
      c.y[i]=y;
      c.rx[i]=rx;
      c.ry[i]=ry;
      c.ra[i]=s.ra;
      c.de[i]=s.de;
      c.mag[i]=s.mag;
      c.select[i]=0;

      if (i>=NMAX)
	break;
      i++;
    }
  }
  fclose(file);
  c.n=i;

  return c;
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
    forward(img.ra0,img.de0,s.ra,s.de,&rx,&ry);
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
