#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cel.h"
#include "cpgplot.h"
#include "qfits.h"
#include <gsl/gsl_multifit.h>
#include <getopt.h>

#define LIM 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define NMAX 4096
#define MAXPOINTERR 40        // Maximum pointing error (pixels)
                              // Do not auto-choose reference stars closer than this to any edge (in pixels) because matching star may
                               // fall outside astrometric field. This should change to arcseconds in the future and be used to widen
                               // the astrometric catalog so that it's scale is grater than that of the imaged catalog just enough
                               // to be sure it includes all imaged stars even at the worst pointing error.
#define MAXROT  10*D2R        // Maximum expected rotation error (radians)
#define DISTMATCHTOL 4
#define MAGMATCHTOL 1
#define MAXSCALERR 1.05
#define MAXMAGERR 4

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
struct map {
  double lat,lng;
  float alt;
  int site_id;
  char observer[32];
} m;

struct image read_fits(char *filename,int pnum);
int fgetline(FILE *,char *,int);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double,double,double,double,double *,double *);
void lfit2d(float *x,float *y,float *z,int n,float *a);
struct catalog read_pixel_catalog(char *filename);
double gmst(double mjd);
double modulo(double x,double y);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);
double sex2dec(char *s);

// Read astrometric catalog
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
    fprintf(stdout,"%s not found!\n",filename);
    exit(0);
  }
  while (!feof(file)) {
    fread(&s,sizeof(struct star),1,file);
    if (s.mag>mmin)
      continue;
    precess(mjd0,s.ra,s.de,img.mjd,&ra,&de);
    forward(img.ra0,img.de0,ra,de,&rx,&ry);
    x=img.x0+1.0/sx*(cos(angle*D2R)*rx+sin(angle*D2R)*ry);
    y=img.y0+1.0/sy*(-sin(angle*D2R)*rx+cos(angle*D2R)*ry);
      /*
    } else if (t.state==1) {
      dx=rx-t.a[0];
      dy=ry-t.b[0];
      d=t.a[1]*t.b[2]-t.a[2]*t.b[1];
      x=(t.b[2]*dx-t.a[2]*dy)/d;
      y=(t.a[1]*dy-t.b[1]*dx)/d;
    }
    */
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
    precess(mjd0,s.ra,s.de,img.mjd,&ra,&de);
    forward(img.ra0,img.de0,ra,de,&rx,&ry);
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

// Get observing site
void get_site(int site_id)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64],filename[LIM],*env;

  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/sites.txt",env);
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with site information not found!\n");
    return;
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip
    if (strstr(line,"#")!=NULL)
      continue;

    // Strip newline
    line[strlen(line)-1]='\0';

    // Read data
    sscanf(line,"%4d %2s %lf %lf %f",
	   &id,abbrev,&lat,&lng,&alt);
    strcpy(observer,line+38);

    // Change to km
    alt/=1000.0;
    
    if (id==site_id) {
      m.lat=lat;
      m.lng=lng;
      m.alt=alt;
      m.site_id=id;
      strcpy(m.observer,observer);
    }

  }
  fclose(file);

  return;
}

// Identify matching triangles in imaged and astrometric catalogs 
void identify_triangles(struct catalog *cat, struct catalog *ast, int *nselect, int tricat[3], int triast[3])
{
  float r,d;
  float mmax,mmin,matchscale;
  float error=MAXPOINTERR;   // maximum expected pointing error in pixels
  int i,s1,s2,s3,m1,m2,m3;
  static int nmatch;
  static float mave,size;
  
    
  // Form a triangle of reference stars from imaged catalog. The following criteria is used:
  // Triangle should cover as much area of FOV as possible.
  // Cannot choose stars too near edges because they may fall outside of astrometric catalog
  // 1st star will be chosen with lower than average magnitude and separated from any edge as
  // much as estimated pointing error. Estimated error is given in pixels for now.
  // 2nd and 3rd stars will also be chosen with lower than average magnitude.
  //
  // Acceptable minimum distance from stars in triangle is set as image width minus pointing error
  // In the future min distance between stars and estimated pointing error should be given in arcsecs and converted
  // to pixels
  
  // if no reference triangle has been previousely selected begin star selection from max size and min magnitude
  if(*nselect==0){
    size=img.naxis2-4*error;
    // Compute magnitude average
    mave=0;
    for (i=0;i<cat->n;i++) {
      r=cat->mag[i];
   	  if(r > mmax)
   	    mmax=r;
   	  if(r < mmin)
   	    mmin=r;
   	  mave+=r;
    }
    mave/=cat->n;
  }

  if(size<0)
    fprintf(stdout,"Error: Pointing inaccuracy too high\n");
  
  fprintf(stdout,"Stars in imaged catalog: %d, stars in astrometric catalog: %d.\n", cat->n, ast->n);
  fprintf(stdout,"MMin:%.1f,MMax:%.1f,MAve.:%.1f \n",mmin,mmax,mave);

  // if last time a suitable reference triangle was found
  if(*nselect==3){
    // and also a match was found
    if(nmatch==3){
      // this time search for another match beginning with the following 3rd star
      triast[2]++;
      fprintf(stdout,"Try different match.\n");      
    }
    else{
      // no match was found for last suitable triangle
      // so search nor another reference triangle starting from the following 3rd  star
      tricat[2]++;
      // erase matching triangle
      triast[0]=0;
      triast[1]=0;
      triast[2]=0;
      fprintf(stdout,"Try different triangle.\n");      
    }
  }
  
  *nselect=0;  
  // search for candidate triangles satisfying size minimum and stars magnitude, later we can fall back to smaller size if no suitable triangle
  // is found
  while((*nselect<3) && (size > 10)){
  
    fprintf(stdout,"Looking for triangles with min size %4.0f\n",size);
    
    // Search candidate for 1st star
    
    for(s1=tricat[0];(*nselect<3) && (s1 < cat->n);s1++){
      // discard if outside high magnitude range
      r=cat->mag[s1];
      if(r > mave){
        fprintf(stdout,"Star %d discarded for high magnitude.\n",s1);
        continue;
      }
      // discard if too close to any edge
      if((cat->x[s1] < error) || (cat->x[s1] > (img.naxis1-error))){
        fprintf(stdout,"Star %d discarded for close to edge.\n",s1);
        continue;
      }
      if((cat->y[s1] < error) || (cat->y[s1] > (img.naxis2-error))){
        fprintf(stdout,"Star %d discarded for close to edge.\n",s1);
        continue;
      }
      // s1 is viable candidate for 1st star
      fprintf(stdout,"Candidate star 1: %d, magnitude: %.1f \n",s1,r);
      tricat[0]=s1;
      *nselect=1;
      
      fprintf(stdout,"Looking for 2nd star of triangle.\n");
      // Search candidate for 2nd star
      for(s2=tricat[1];(*nselect<3) && (s2 < cat->n);s2++){
        // discard candidate if outside magnitude range
        r=cat->mag[s2];
        if(r > mave){
          fprintf(stdout,"Star %d discarded for high magnitude.\n",s2);
          continue;
        }
        // discard if too close to any edge
        if((cat->x[s2] < error) || (cat->x[s2] > (img.naxis1-error))){
          fprintf(stdout,"Star %d discarded for close to edge.\n",s2);
          continue;
        }
        if((cat->y[s2] < error) || (cat->y[s2] > (img.naxis2-error))){
          fprintf(stdout,"Star %d discarded for close to edge.\n",s2);
          continue;
        }
        // discard if too close to star 1
        r=sqrt(pow(cat->x[s2]-cat->x[s1],2)+pow(cat->y[s2]-cat->y[s1],2));
        if(r < size){
          fprintf(stdout,"Star %d discarded for close to 1st star.\n",s2);
          continue;
        }
        // s2 is viable candidate for 2nd star
        fprintf(stdout,"Candidate star 2: %d, magnitude: %.1f, distance to star 1: %.1f \n",s2,cat->mag[s2],r);
        tricat[1]=s2;
        *nselect=2;
        
        fprintf(stdout,"Looking for 3rd star of triangle.\n");        
        // Search candidate for 3rd star beginning from the last candidate found
        for(s3=tricat[2];(*nselect<3) && (s3 < cat->n);s3++){
          // discard candidate if outside magnitude range
          r=cat->mag[s3];
          if(r > mave){
            fprintf(stdout,"Star %d discarded for high magnitude.\n",s3);
            continue;
          }
          // discard if too close to any edge
          if((cat->x[s3] < error) || (cat->x[s3] > (img.naxis1-error))){
            fprintf(stdout,"Star %d discarded for close to edge.\n",s3);
            continue;
          }
          if((cat->y[s3] < error) || (cat->y[s3] > (img.naxis2-error))){
            fprintf(stdout,"Star %d discarded for close to edge.\n",s3);
            continue;
          }
          // discard if too close to star 1
          r=sqrt(pow(cat->x[s3]-cat->x[s1],2)+pow(cat->y[s3]-cat->y[s1],2));
          if(r < size){
            fprintf(stdout,"Star %d discarded for close to 1st star.\n",s3);
            continue;
          }
          // discard if too close to star 2
          r=sqrt(pow(cat->x[s3]-cat->x[s2],2)+pow(cat->y[s3]-cat->y[s2],2));
          if(r < size){
            fprintf(stdout,"Star %d discarded for close to 2nd star.\n",s3);
            continue;
          }
          // s3 is viable candidate for 3rd star
          fprintf(stdout,"Candidate star 3: %d, magnitude: %.1f, distance to star 2: %.1f \n",s3,cat->mag[s3],r);
          tricat[2]=s3;
          *nselect=3;
          
          
          // ****************************
          // ****************************
          // ****************************
          // ****************************
          //Ref.triangle: 8, 86, 91, 
          //Match triangle: 125, 64, 115,                    
          //Ref.triangle: 8, 94, 161,
          //Match triangle: 125, 1, 90,
//Ref.triangle: 8, 159, 161,
//Match triangle: 125, 55, 90,
//Ref.triangle: 8, 152, 129,
//Match triangle: 151, 54, 120,
//Ref.triangle: 8, 77, 129,
//Match triangle: 125, 24, 120,

//Está matcheando!!!
//ahora habría que descartar mirrors verticales y/o horizontales y limitar la rotación para llegar más directo

          if(s1==8 && s2==86 && s3==91){
            fprintf(stdout,"Special reference triangle found!!!\n");
          }
          
          // ************************************************************
          // Calibration triangle candidate found!
          // time to look for matching triangle in astrometric catalog...
          // ************************************************************
                 
          // at this point we got 3 reference stars forming suitable triangle, will look for matching triangle in astrometric catalog
          // will search by testing relative magnitude and relative distance between suspect matching stars
          // a scale tolerance is given as a ratio for distances (triangle size) and is assumed as affecting the whole image,
          // this tolerance could be high.
          // Also a magnitude tolerance is defined as an absolute value, this should not be too high (1..2 magnitudes)
          // 
          // relative distance between stars of astrometric catalog should closely match relative distances between stars of imaged catalog
          // a tolerance for this match is given also as a ratio and should be low
          // also magnitude difference between stars of astrometric cat should closely match magn diff between stars of img cat.
          // a tolerance for this magnitude match is given and should not be high (0.5 magn)
          
          nmatch=0;
          // Search for candidate 1st match star
          for(m1=triast[0];((nmatch<3) && (m1<ast->n));m1++){
            // discard if outside magnitude tolerance
            // this tolerance includes magnitude extraction tolerance on sextractor tool 
            r=fabs(cat->mag[s1] - ast->mag[m1]);
            if(r > MAXMAGERR){
              fprintf(stdout,"Match 1st star %d discarded for mag. difference.\n",m1);
              continue;
            }
            r=ast->mag[m1];
            // m1 is viable candidate for match to 1st star
            fprintf(stdout,"Candidate match star 1: %d, magnitude: %.1f \n",m1,r);
            triast[0]=m1;
            nmatch=1;
            // Search candidate for 2nd star
            for(m2=triast[1];(nmatch<3) && (m2<ast->n);m2++){
              // magn difference between star 1 and star 2 should be very close to magn diff between matching star 1 and matching star 2
              // discard if difference is outside match tolerance
              r=fabs(cat->mag[s2] - cat->mag[s1]);
              if((r - fabs(ast->mag[m2] - ast->mag[m1])) > MAGMATCHTOL){
                fprintf(stdout,"Match 2nd star %d discarded for mag. match tolerance.\n",m2);
                continue;
              }
              // distance from star 1 to star 2 should resemble distance from matching star 1 to matching star 2 allowing for
              // scale tolerance
              r=sqrt(pow(cat->x[s2]-cat->x[s1],2)+pow(cat->y[s2]-cat->y[s1],2));
              // discard if distance between m2 and m1 is too far from distance between s2 and s1
              if(sqrt(pow(ast->x[m2]-ast->x[m1],2)+pow(ast->y[m2]-ast->y[m1],2)) / r > MAXSCALERR){
                fprintf(stdout,"Match 2nd star %d discarded for distance to 1st match star error.\n",m2);
                continue;
              }
              if(r / sqrt(pow(ast->x[m2]-ast->x[m1],2)+pow(ast->y[m2]-ast->y[m1],2)) > MAXSCALERR){
                fprintf(stdout,"Match 2nd star %d discarded for distance to 1st match star error.\n",m2);
                continue;
              }
              // check angle
              d=atan((cat->y[s2]-cat->y[s1]) / (cat->x[s2]-cat->x[s1]));
              if(fabs(d - atan((ast->y[m2]-ast->y[m1]) / (ast->x[m2]-ast->x[m1]))) > MAXROT){
                fprintf(stdout,"Match 2nd star %d discarded for rotation error.\n",m2);
                continue;
              }                                          
              // keep scale to check for coherence with third matching star
              matchscale=sqrt(pow(ast->x[m2]-ast->x[m1],2)+pow(ast->y[m2]-ast->y[m1],2)) / r;
              fprintf(stdout,"Triangle size match scale: %f.2.\n",matchscale);
              // m2 is viable candidate for match to 2nd star
              r=ast->mag[m2];
              fprintf(stdout,"Candidate match star 2: %d, magnitude: %.1f \n",m2,r);

              fprintf(stdout,"Angle star 1 to star 2: %f \n",d);
              d=atan((ast->y[m2]-ast->y[m1]) / (ast->x[m2]-ast->x[m1]));
              fprintf(stdout,"Angle match 1 to match 2: %f \n",d);
              d=fabs(d - atan((ast->y[m2]-ast->y[m1]) / (ast->x[m2]-ast->x[m1])));
              fprintf(stdout,"Angle error %f \n",d);
              
              triast[1]=m2;
              nmatch=2;
              // Search candidate for 2nd star
              // start from the following star that was found as viable candidate last time
              // this allows for selecting a new match each time operator launches auto-calibrate
              for(m3=triast[2];(nmatch<3) && (m3<ast->n);m3++){
                // discard if magn difference is outside match tolerance
                r=fabs(cat->mag[s3] - cat->mag[s1]);
                if((r - fabs(ast->mag[m3] - ast->mag[m1])) > MAGMATCHTOL){
                  fprintf(stdout,"Match 3rdd star %d discarded for mag. match tolerance.\n",m3);
                  continue;
                }
                // scale was fixed when found 2nd match candidate, no different scale error is allowed for 3rd matching star 
                // discard if distance between m3 and m1 is not very similar to distance between s3 and s1
                //**************se esta descartando el triangulo que corresponde por distancia, aumentar tolerancia y revisar...
                r=matchscale * sqrt(pow(cat->x[s3]-cat->x[s1],2)+pow(cat->y[s3]-cat->y[s1],2));
                d=sqrt(pow(ast->x[m3]-ast->x[m1],2)+pow(ast->y[m3]-ast->y[m1],2));                
                if(fabs(d - r) > DISTMATCHTOL){
                  fprintf(stdout,"Match 3rd star %d discarded for distance to 1st match star error: %f against %f.\n",m3,r,d);
                  continue;
                }
                // discard if distance between m3 and m2 is not very similar to distance between s3 and s2
                r=matchscale * sqrt(pow(cat->x[s3]-cat->x[s2],2)+pow(cat->y[s3]-cat->y[s2],2));
                if(fabs(sqrt(pow(ast->x[m3]-ast->x[m2],2)+pow(ast->y[m3]-ast->y[m2],2)) - r) > DISTMATCHTOL){
                  fprintf(stdout,"Match 3rd star %d discarded for distance to 2nd match star error.\n",m3);
                  continue;
                }
                // check angle
                d=atan((cat->y[s3]-cat->y[s1]) / (cat->x[s3]-cat->x[s1]));
                if(fabs(d - atan((ast->y[m3]-ast->y[m1]) / (ast->x[m3]-ast->x[m1]))) > MAXROT){
                  fprintf(stdout,"Match 3rd star %d discarded for rotation error.\n",m3);
                  continue;
                }
                // m3 is viable candidate for match to 3rd star
                r=ast->mag[m3];
                fprintf(stdout,"Candidate match star 3: %d, magnitude: %.1f \n",m3,r);
                triast[2]=m3;
                nmatch=3;
              }
              // if no match found next try will search again for all possible 3rd match stars
              if(nmatch<3) triast[2]=0;
            }
            // if no match found next try will search again for all possible 2nd match stars
            if(nmatch<3) triast[1]=0;
          }
          if(nmatch<3){
            triast[0]=0;
            triast[1]=0;
            triast[2]=0;
            // try new 3rd star
            *nselect=2;
            fprintf(stdout,"No match found for this triangle\n");
          }
        }
        // if no triangle selected we will try next 2nd star and try again all possible 3rd stars
        if(*nselect<3) tricat[2]=0;
        fprintf(stdout,"Try next 2nd star.\n");
      }
      // if no triangle selected we will try next 1st star and try again all possible 2nd stars
      if(*nselect<3) tricat[1]=0;
      fprintf(stdout,"Try next 1st star.\n");
    }
    // at this point either found 3 candidate stars or already tried all candidate stars without success
    // if loop continues will look for smaller triangles and accept higher magnitudes next iteration
    if(*nselect<3){
      size/=1.25;
      mave+=0.25;    
    }
    fprintf(stdout,"Decrase minimum size and increase max magnitude for next search.\n");
  }
  // if we could not find suitable triplet erase tricat[] 
  if(*nselect<3){
    tricat[0]=0;
    tricat[1]=0;
    tricat[2]=0;
    fprintf(stdout,"Did not find any suitable reference triangle.\n");
  }
  return;
}


int main(int argc,char *argv[])
{
  int i;
  float xmin,xmax,ymin,ymax,zmin,zmax,width;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float x,y,r,rmin=1.0,rmax=10.0,mmin=5.0,mmax=10.0;
  struct catalog cat,ast;
  char c;
  int redraw=1,click=0,nselect=0,plotstars=1;
  char filename[128],sra[20],sde[20],cam[15],mount[15];
  float h,q,sx,sy,mag=9,fw,fh;
  int nx,ny;
  FILE *file;
  char *env,starfile[128];
  int tricat[3],triast[3],n;

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(starfile,"%s/data/tycho2.dat",env);

  // Geographic position
  env=getenv("ST_COSPAR");
  get_site(atoi(env));

  // Read image
  img=read_fits(argv[1],0);
  sprintf(filename,"%s.cat",argv[1]);

  printf("Image read\n");

  // Initial transformation
  if (argc==7) {
    sx=-atof(argv[2]);
    sy=-sx;
    img.ra0=atof(argv[3]);
    img.de0=atof(argv[4]);
    q=atof(argv[5]);
    mag=atof(argv[6]);
  } else {
    file=fopen("position.txt","r");
    if (file==NULL) {
      fprintf(stdout,"No position file found\n");
      return 0;
    }
    fscanf(file,"%s %s",sra,sde);
    fclose(file);
    
    // Get parameters
    img.ra0=15.0*sex2dec(sra);
    img.de0=sex2dec(sde);

    // Hour angle and parallactic angle
    h=gmst(img.mjd)+m.lng-img.ra0;
    q=atan2(sin(h*D2R),(tan(m.lat*D2R)*cos(img.de0*D2R)-sin(img.de0*D2R)*cos(h*D2R)))*R2D;
    printf("Hour angle: %.3f deg, parallactic angle: %.3f deg\n",h,q);

    // Get pixel scale params from camera file
    file=fopen("camera.txt","r");
    if (file==NULL) {
      sx=-36.15;
      sy=33.22;
      fprintf(stdout,"No camera file found, using default pixel scale values %3.2f %3.2f\n",sx,sy);
    }
    else{
      // Obtain FOV and image resolution from camera file
      fscanf(file,"%s %f %f %d %d %s",cam,&fw,&fh,&nx,&ny,mount);
      fclose(file);
      sx=fw/nx*3600;
      sy=fh/ny*3600;
      // Check scheduled camera resolution against FITS image file resolution
      if((abs(nx)!=img.naxis1) || (abs(ny)!=img.naxis2)){
        fprintf(stdout,"Warning: scheduled camera resolution %dx%d does not match image resolution %dx%d\n",nx,ny,img.naxis1,img.naxis2);
      }
    }

  }
  img.x0=0.5*(float) img.naxis1;
  img.y0=0.5*(float) img.naxis2;

  // Read catalogs
  cat=read_pixel_catalog(filename);
  ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);

  // Open PGPlot server
  cpgopen("/xs");
//  cpgwnad(0.0,img.naxis1,0.0,img.naxis2);
//  cpgsfs(2);
//  cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
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


  // For ever loop
  for (;;) {      
    if (redraw==1) {
      cpgeras();

      cpgsvp(0.1,0.95,0.1,0.95);
      cpgwnad(xmin,xmax,ymin,ymax);
      cpglab("x (pix)","y (pix)"," ");
      cpgsfs(2);
      cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      cpgimag(img.z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,img.zmin,img.zmax,tr);
      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

      
      // Plot triangle of calibration stars
      if(nselect>=2){
        // imaged catalog
        cpgmove(cat.x[tricat[0]],cat.y[tricat[0]]);
        fprintf(stdout,"\nRef.triangle: ");
        fprintf(stdout,"%d, ",tricat[0]);
        for(i=1;(i<3 && i<nselect);i++){
          cpgdraw(cat.x[tricat[i]],cat.y[tricat[i]]);
          fprintf(stdout,"%d, ",tricat[i]);
        }
        cpgdraw(cat.x[tricat[0]],cat.y[tricat[0]]);
        // astrometric catalog
        if(triast[2]>0){
          cpgmove(ast.x[triast[0]],ast.y[triast[0]]);
          fprintf(stdout,"\nMatch triangle: ");
          fprintf(stdout,"%d, ",triast[0]);
          for(i=1;(i<3 && i<nselect);i++){
            cpgdraw(ast.x[triast[i]],ast.y[triast[i]]);
            fprintf(stdout,"%d, ",triast[i]);
          }
          cpgdraw(ast.x[triast[0]],ast.y[triast[0]]);
        }
      }
      fflush(stdout);
      
      // Plot catalogs
      if (plotstars==1) {
	      cpgsci(3);
	      for (i=0;i<cat.n;i++) {
         	  r=rmax-(rmax-rmin)*(cat.mag[i]-mmin)/(mmax-mmin);
	        r*=img.naxis1/752.0;
	        //r*=0.1;
	        if ((n=cat.select[i])!=0){
	          cpgpt1(cat.x[i],cat.y[i],0);
                }
	        //else
	          //cpgpt1(cat.x[i],cat.y[i],4);
	        cpgcirc(cat.x[i],cat.y[i],r);
	      }
      }
      
      cpgsci(4);
      for (i=0;i<ast.n;i++) {
	      r=rmax-(rmax-rmin)*(ast.mag[i]-mmin)/(mmax-mmin);
	      // Upscale for image size
	      r*=img.naxis1/752.0;
	      if (ast.select[i]!=0)
	        cpgpt1(ast.x[i],ast.y[i],0);
        cpgcirc(ast.x[i],ast.y[i],r);
      }
      cpgsci(1);
      redraw=0;
    }

    cpgcurs(&x,&y,&c);

    // Quit
    if (c=='q')
      break;

    // Fit
    if (c=='f' && nselect>=3) {
      fit_transformation(cat,ast,nselect);
      ast=reread_astrometric_catalog(starfile,mag+1);
      tricat[0]=0;
      triast[0]=0;
      tricat[1]=0;
      triast[1]=0;
      tricat[2]=0;
      triast[2]=-1;
      nselect=0;
      redraw=1;
    }

    // Reread
    if (c=='R') {
      for (i=0;i<cat.n;i++) {
        cat.select[i]=0;
      }
      for (i=0;i<ast.n;i++) {
        ast.select[i]=0;
      }
      for(i=0;i<3;i++){
        tricat[i]=0;
        triast[i]=0;
        img.a[i]=0;
        img.b[i]=0;
      }
      triast[2]=-1;
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
      //ast=reread_astrometric_catalog(starfile,mag);
      redraw=1;
      nselect=0;
    }

    // Select pixel catalog
    if (c=='a' && click==0) {
      i=select_nearest(cat,x,y);
      fprintf(stdout,"%d",i);      
      cat.select[i]=nselect+1;
      if(nselect<4) 
        tricat[nselect]=i;
      redraw=1;
      click=1;
    }

    // Select catalog
    if (c=='b' && click==1) {
      i=select_nearest(ast,x,y);
      fprintf(stdout,"%d",i);            
      ast.select[i]=nselect+1;
      if(nselect<4)
        triast[nselect]=i;
      redraw=1;
      click=0;
      nselect++;
    }

    // Select pixel catalog
    if (c=='1') {
      i=select_nearest(cat,x,y);
      fprintf(stdout,"%d",i);      
    }

    // Select catalog
    if (c=='2') {
      i=select_nearest(ast,x,y);
      fprintf(stdout,"%d",i);            
    }

    // Autoidentify triangles in imaged and astrometric catalogs
    if (c=='i') {
      identify_triangles(&cat,&ast,&nselect,tricat,triast);
      redraw=1;
      click=1;
    }

    
    // Plot identified stars
    if (c=='p') {
      if (plotstars==1)
        plotstars=0;
      else if (plotstars==0)
        plotstars=1;
      redraw=1;
    }

    // Match catalogs
    if (c=='m') {
      match_catalogs(&cat,&ast,10.0);
      redraw=1;
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

    // Help
    if (c=='h') {
      printf("Calibrates astrometry. Initially requires manual matching of at least three stars. Use 'a' to select star on the image, then 'b' to select star from the catalog, then 'f' to fit");
      printf("q     Quit\n");
      printf("a     Select star on image\n");
      printf("b     Select star from catalog\n");
      printf("i     Autoselect calibration stars from catalog\n");
      printf("c     Center image on pixel\n");
      printf("f     Fit calibration\n");
      printf("m     Match stars using current calibration\n");
      printf("z/+   Zoom in on cursor\n");
      printf("x/-   Zoom out on cursor\n");
      printf("p     Plot sextractor catalog\n");
      printf("r     Reset zoom\n");
      printf("R     Reset fit and star selections\n");      
    }
    redraw=1;
    continue;
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
  img.mjd=atof(qfits_query_hdr(filename,"MJD-OBS"));

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
  avg=s1/(float) (img.naxis1*img.naxis2);
  std=sqrt(s2/(float) (img.naxis1*img.naxis2)-avg*avg);
  printf("%f %f\n",avg,std);
  img.zmin=avg-4.0*std;
  img.zmax=avg+6.0*std;
  
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

// Get a x and y from a RA and Decl
void forward(double ra0,double de0,double ra,double de,double *x,double *y)
{
  int i;
  char pcode[4]="STG";
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
  *x *=3600.;
  *y *=3600.;

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

// Get a RA and Decl from x and y
void reverse(double ra0,double de0,double x,double y,double *ra,double *de)
{
  int i;
  char pcode[4]="STG";
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

// Read pixel catalog
struct catalog read_pixel_catalog(char *filename)
{
  int i=0;
  FILE *file;
  char line[LIM];
  struct catalog c;

  // Read catalog
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stdout,"%s not found!\n",filename);
    exit(0);
  }
  while (fgetline(file,line,LIM)>0) {
    if (strstr(line,"#")!=NULL) 
      continue;
    sscanf(line,"%f %f %f",&c.x[i],&c.y[i],&c.mag[i]);
    c.select[i]=0;
    i++;
  }
  fclose(file);
  c.n=i;

  return c;
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
