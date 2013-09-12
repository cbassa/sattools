#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qfits.h"
#include "cpgplot.h"
#include "cel.h"
#include <gsl/gsl_multifit.h>

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
  float avg,std;
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
int fgetline(FILE *file,char *s,int lim);
int select_nearest(struct catalog c,float x,float y);

void plot_objects(char *filename)
{
  int i;
  FILE *file;
  float x0,y0,x1,y1,texp;
  int id;
  char line[LIM],catalog[128],dummy[128],text[8];

  file=fopen(filename,"r");
  while (fgetline(file,line,LIM)>0) {
    sscanf(line,"%s %f %f %f %f %f %d %s",dummy,&x0,&y0,&x1,&y1,&texp,&id,catalog);

    if (strstr(catalog,"classfd")!=NULL)
      cpgsci(4);
    else
      cpgsci(1);

    cpgmove(x0,y0);
    cpgdraw(x1,y1);
    cpgpt1(x0,y0,4);
    sprintf(text," %05d",id);
    cpgtext(x0,y0,text);
  }
  fclose(file);
  cpgsci(1);

  return;
}

int main(int argc,char *argv[])
{
  int i;
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
  float sx,sy,q;
  char *env,idfile[128];
  float r,rmin=1.0,rmax=10.0,mmin=1.0,mmax=8.0,mag=8.0;

  // Environment variables
  env=getenv("ST_DATADIR");

  // Read image
  img=read_fits(argv[1],0);

  sprintf(idfile,"%s.id",argv[1]);

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
      cpgsci(1);
      
      cpgsvp(0.1,0.95,0.1,0.95);
      cpgwnad(xmin,xmax,ymin,ymax);
      cpglab("x (pix)","y (pix)"," ");
      cpgsfs(2);
      cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);

      cpgimag(img.z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

      plot_objects(idfile);

      redraw=0;
    }

    // Get cursor
    cpgcurs(&x,&y,&c);

    // Quit
    if (c=='q')
      break;


    // Center
    if (c=='c') {
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // Zoom
    if (c=='z' || c=='+') {
      width/=1.5;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // Unzoom
    if (c=='x' || c=='+' || c=='=') {
      width*=1.5;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
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

