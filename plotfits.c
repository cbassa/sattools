#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cel.h"
#include "cpgplot.h"
#include "qfits.h"
#include <gsl/gsl_multifit.h>
#include <getopt.h>
#include <time.h>

#define LIM 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define NMAX 4096
#define MAXPOINTERR 200        // Maximum pointing error (pixels)
                              // Algorithm will not choose reference stars closer than this amount to any edge (in pixels)
                              // because matching star may fall outside astrometric field.
                              // This should change to arcseconds in the future and be used to widen
                               // the astrometric catalog so that it's scale is grater than that of the imaged catalog just enough
                               // to be sure it includes all imaged stars even at the worst pointing error.
#define MAXROT  10*D2R        // Maximum expected rotation error (radians)
#define MAXSCALERR 1.05       // Expected image to astrometric map scaling error
#define DISTMATCHTOL 6        // Distance tolerance in pixels between matching stars after aplying scale and rotation
#define MAXMAGERR 3         // Expected magnitude error between imaged stars and corresponding astometric catalog stars
#define MAGMATCHTOL 1       // Relative magnitude between matching stars tolerance
#define DEFMATCHVALRATIO 0.3     // Default ratio of imaged stars that must fit into astrometric catalog after applying matching transformation (can be adjusted at runtime)
#define	AUTOMAGLIM 1				// Automatically set astrometric catalog magnitude limit
#define DEBUG 1

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
float matchvalratio=DEFMATCHVALRATIO;


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
  float r,d,d2,x,y;
  float mag12,mag13,mag23,ang12,ang13,dis12,dis13,dis23,cenX,cenY;
  float mmax,mmin,matchscale,matchrot,matchXtras,matchYtras;
  float error=MAXPOINTERR;   // maximum expected pointing error in pixels
	int i,j,n,n2,s1,s2,s3,m1,m2,m3;
  static int nmatch;
  clock_t start_t, end_t;
	static float mave,size;
  
  start_t=clock();
    
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

  
  // if a reference triangle match was already found, we want to find another one, so we don't start
  // from zero
  if((triast[0]!=0) || (triast[1]!=0) || (triast[2]!=0)){
    // this time search for another match beginning with the following 3rd star
    triast[2]++;
#if DEBUG>0
    fprintf(stdout,"Try different match.\n");
#endif
  }
  else{
    // Either no match was found or this is the first time the routine is launched so
    // we calculate reference triangle parameters...
	//  size=(img.naxis2-error)/2;
		size=(img.naxis2-error)/3;
		// Compute magnitude average
		mave=0;
		mmax=0;
		mmin=10;
		for (i=0;i<cat->n;i++) {
			r=cat->mag[i];
	 	  if(r > mmax)
	 	    mmax=r;
	 	  if(r < mmin)
	 	    mmin=r;
	 	  mave+=r;
		}
		mave/=cat->n;
	#if DEBUG>2
		fprintf(stdout,"Image height: %d, Expected error: %f, Ref.Tria.Size: %f\n",img.naxis2,error,size);      
	#endif

		if(size<0){
			fprintf(stdout,"Warning: Pointing inaccuracy too high\n");
			size=10;
		}
	#if DEBUG>0
		fprintf(stdout,"Stars in imaged catalog: %d, stars in astrometric catalog: %d.\n", cat->n, ast->n);
		fprintf(stdout,"MMin:%2.1f,MMax:%2.1f,MAve.:%2.1f \n",mmin,mmax,mave);
	#endif
    
    triast[0]=0;
    triast[1]=0;
    triast[2]=0;
#if DEBUG>2
    fprintf(stdout,"Try different triangle. Size: %f\n",size);      
#endif
  }
 
 
  *nselect=0;  
  nmatch=0;
  // search for candidate triangles satisfying size minimum and stars magnitude, later we can fall back to smaller size if no suitable triangle
  // is found
//  while((*nselect<3) && (size > 10)){
  while((nmatch<3) && (size > 10)){
//	  *nselect=0;  
  
#if DEBUG>1
    fprintf(stdout,"Looking for triangles with min size %4.0f\n",size);
#endif    
    // Search candidate for 1st star
    
    for(s1=tricat[0];(nmatch<3) && (s1 <= cat->n);s1++){
      // discard if outside high magnitude range
      r=cat->mag[s1];
      if(r > mave){
#if DEBUG>2
        fprintf(stdout,"Star %d discarded for high magnitude.\n",s1);
#endif
        continue;
      }
      // discard if too close to any edge
      if((cat->x[s1] < error/2) || (cat->x[s1] > (img.naxis1-error/2))){
#if DEBUG>2
        fprintf(stdout,"Star %d discarded for close to edge.\n",s1);
#endif
        continue;
      }
      if((cat->y[s1] < error/2) || (cat->y[s1] > (img.naxis2-error/2))){
#if DEBUG>2
        fprintf(stdout,"Star %d discarded for close to edge.\n",s1);
#endif
        continue;
      }
      // s1 is viable candidate for 1st star
#if DEBUG>1
      fprintf(stdout,"Candidate star 1: %d, magnitude: %.1f \n",s1,r);
#endif
      tricat[0]=s1;
      *nselect=1;
      
#if DEBUG>1
      fprintf(stdout,"Looking for 2nd star of triangle.\n");
#endif
      // Search candidate for 2nd star
      for(s2=tricat[1];(nmatch<3) && (s2 <= cat->n);s2++){
        // discard candidate if outside magnitude range
        r=cat->mag[s2];
        if(r > mave){
#if DEBUG>2
          fprintf(stdout,"Star %d discarded for high magnitude.\n",s2);
#endif
          continue;
        }
        // discard if too close to any edge
        if((cat->x[s2] < error/2) || (cat->x[s2] > (img.naxis1-error/2))){
#if DEBUG>2
          fprintf(stdout,"Star %d discarded for close to edge.\n",s2);
#endif
          continue;
        }
        if((cat->y[s2] < error/2) || (cat->y[s2] > (img.naxis2-error/2))){
#if DEBUG>2
          fprintf(stdout,"Star %d discarded for close to edge.\n",s2);
#endif
          continue;
        }
        // discard if too close to star 1
        r=sqrt(pow(cat->x[s2]-cat->x[s1],2)+pow(cat->y[s2]-cat->y[s1],2));
        if(r < size){
#if DEBUG>2
          fprintf(stdout,"Star %d discarded for close to 1st star.\n",s2);
#endif
          continue;
        }
        // s2 is viable candidate for 2nd star
#if DEBUG>1
        fprintf(stdout,"Candidate star 2: %d, magnitude: %.1f, distance to star 1: %.1f \n",s2,cat->mag[s2],r);
#endif
        tricat[1]=s2;
        *nselect=2;
        
#if DEBUG>1
        fprintf(stdout,"Looking for 3rd star of triangle.\n");        
#endif
        // Search candidate for 3rd star beginning from the last candidate found
        for(s3=tricat[2];(nmatch<3) && (s3 <= cat->n);s3++){
          // discard candidate if outside magnitude range
          r=cat->mag[s3];
          if(r > mave){
#if DEBUG>2
            fprintf(stdout,"Star %d discarded for high magnitude.\n",s3);
#endif
            continue;
          }
          // discard if too close to any edge
          if((cat->x[s3] < error/2) || (cat->x[s3] > (img.naxis1-error/2))){
#if DEBUG>2
            fprintf(stdout,"Star %d discarded for close to edge.\n",s3);
#endif
            continue;
          }
          if((cat->y[s3] < error/2) || (cat->y[s3] > (img.naxis2-error/2))){
#if DEBUG>2
            fprintf(stdout,"Star %d discarded for close to edge.\n",s3);
#endif
            continue;
          }
          // discard if too close to star 1
          r=sqrt(pow(cat->x[s3]-cat->x[s1],2)+pow(cat->y[s3]-cat->y[s1],2));
          if(r < size){
#if DEBUG>2
            fprintf(stdout,"Star %d discarded for close to 1st star.\n",s3);
#endif
            continue;
          }
          // discard if too close to star 2
          r=sqrt(pow(cat->x[s3]-cat->x[s2],2)+pow(cat->y[s3]-cat->y[s2],2));
          if(r < size){
#if DEBUG>2
            fprintf(stdout,"Star %d discarded for close to 2nd star.\n",s3);
#endif
            continue;
          }
          // s3 is viable candidate for 3rd star
#if DEBUG>1
          fprintf(stdout,"Candidate star 3: %d, magnitude: %.1f, distance to star 2: %.1f \n",s3,cat->mag[s3],r);
#endif

#if DEBUG>1
          fprintf(stdout,"\n");
          fprintf(stdout,"Reference triangle: %d, %d, %d\n",s1,s2,s3);
          fprintf(stdout,"Magnitudes: %f, %f, %f\n",cat->mag[s1],cat->mag[s2],cat->mag[s3]);          
#endif
          tricat[2]=s3;
          *nselect=3;
                    
          // ************************************************************
          // Calibration triangle candidate found!
          // now look for matching triangle in astrometric catalog...
          // ************************************************************
                 
          // at this point we have got 3 reference stars forming suitable triangle, will look for matching triangle in astrometric catalog
          // will search by testing relative magnitude and relative distance between suspect matching stars
          // a scale tolerance is given as a ratio for distances (triangle size) and is assumed as affecting the whole image,
          // this tolerance could be high.
          // Also a magnitude tolerance is defined as an absolute value, this should not be too high (1..2 magnitudes)
          // 
          // relative distance between stars of astrometric catalog should closely match relative distances between stars of imaged catalog
          // a tolerance for this match is given also as a ratio and should be low
          // also magnitude difference between stars of matching triangle should closely match magn diff between stars of imaged
          // reference triangle, tolerance for this magnitude match should not be high (0.5 magn or so)


          nmatch=0;
          // Search for candidate 1st match star
          for(m1=triast[0];((nmatch<3) && (m1<=ast->n));m1++){
            // discard if outside magnitude tolerance
            // this tolerance includes magnitude extraction tolerance on sextractor tool 
            r=fabs(cat->mag[s1] - ast->mag[m1]);
            if(r > MAXMAGERR){
#if DEBUG>2
              fprintf(stdout,"Match 1st star %d discarded for mag. difference.\n",m1);
#endif
              continue;
            }
            r=ast->mag[m1];
            // m1 is viable candidate for match to 1st reference star
#if DEBUG>1
            fprintf(stdout,"Candidate match star 1: %d, magnitude: %.1f \n",m1,r);
#endif
            triast[0]=m1;
            nmatch=1;
            // Calculate some values that will be repeatedly used during match check loop
            mag12=cat->mag[s2] - cat->mag[s1];
            dis12=sqrt(pow(cat->x[s2]-cat->x[s1],2)+pow(cat->y[s2]-cat->y[s1],2));
            ang12=atan((cat->y[s2]-cat->y[s1]) / (cat->x[s2]-cat->x[s1]));
            if((cat->x[s2]-cat->x[s1]) < 0){
              if(ang12 >= 0) ang12 -= M_PI;
              else ang12 += M_PI;
            }
            // Search candidate for 2nd star
            for(m2=triast[1];(nmatch<3) && (m2<=ast->n);m2++){
              // magn difference between ref stars 1 and 2 should be very close to magn diff between matching stars 1 and 2
              // discard if difference is outside match tolerance
              if(fabs(mag12 - (ast->mag[m2] - ast->mag[m1])) > MAGMATCHTOL){
#if DEBUG>2
                fprintf(stdout,"Match 2nd star %d discarded for mag. match tolerance.\n",m2);
#endif
                continue;
              }
              // distance between ref stars 1 and 2 should resemble distance from matching stars 1 and star 2 allowing for
              // scale tolerance
              d=sqrt(pow(ast->x[m2]-ast->x[m1],2)+pow(ast->y[m2]-ast->y[m1],2)) / dis12;
              // discard if distance between m2 and m1 is too far from distance between s2 and s1
              if((d > MAXSCALERR) || (1/d > MAXSCALERR)){
#if DEBUG>2
                fprintf(stdout,"Match 2nd star %d discarded for distance to 1st match star error.\n",m2);
#endif
                continue;
              }
              // keep scale to later check for coherence with third matching star
              // this scale was obtained from distance ratio between ref stars 1 & 2 and match stars 1 & 2
              matchscale=d;
#if DEBUG>1
              fprintf(stdout,"Triangle size match scale: %f.2.\n",matchscale);
#endif
              // check angle
              d2=atan((ast->y[m2]-ast->y[m1]) / (ast->x[m2]-ast->x[m1]));
              if((ast->x[m2]-ast->x[m1]) < 0){
                if(d2 >= 0) d2 -= M_PI;
                else d2 += M_PI;
              }
              if(fabs(d2 - ang12) > MAXROT){
#if DEBUG>2
                fprintf(stdout,"Match 2nd star %d discarded for rotation error.\n",m2);
#endif
                continue;
              }
              // keep rotation for later checks
              // this rotation was obtained from angle between ref stars 1 & 2 and match stars 1 & 2
              matchrot=d2-d;
              // m2 is viable candidate for match to 2nd star
              r=ast->mag[m2];
#if DEBUG>1
              fprintf(stdout,"Candidate match star 2: %d, magnitude: %.1f \n",m2,r);
#endif
#if DEBUG>2
              fprintf(stdout,"Angle star 1 to star 2: %f \n",d);
              d2=atan((ast->y[m2]-ast->y[m1]) / (ast->x[m2]-ast->x[m1]));
              fprintf(stdout,"Angle match 1 to match 2: %f \n",d2);
              d=fabs(d - d2);
              fprintf(stdout,"Angle error %f \n",d);
#endif              
              triast[1]=m2;
              nmatch=2;
              // calculate values that will be used repeatedly
              mag13=cat->mag[s3] - cat->mag[s1];
              mag23=cat->mag[s3] - cat->mag[s2];
              cenX=(cat->x[s1] + cat->x[s2] + cat->x[s3]) / 3;
              cenY=(cat->y[s1] + cat->y[s2] + cat->y[s3]) / 3;
              dis13=matchscale * sqrt(pow(cat->x[s3]-cat->x[s1],2)+pow(cat->y[s3]-cat->y[s1],2));
              dis23=matchscale * sqrt(pow(cat->x[s3]-cat->x[s2],2)+pow(cat->y[s3]-cat->y[s2],2));
              ang13=atan((cat->y[s3]-cat->y[s1]) / (cat->x[s3]-cat->x[s1]));
              if((cat->x[s3]-cat->x[s1]) < 0){
                if(ang13 >= 0) ang13 -= M_PI;
                else ang13 += M_PI;
              }
              // Search candidate for 3rd star
              for(m3=triast[2];(nmatch<3) && (m3<=ast->n);m3++){
                // discard if magn difference is outside match tolerance
                if(fabs(mag13 - (ast->mag[m3] - ast->mag[m1])) > MAGMATCHTOL){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for mag. match tolerance.\n",m3);
#endif
                  continue;
                }
                if(fabs(mag23 - (ast->mag[m3] - ast->mag[m2])) > MAGMATCHTOL){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for mag. match tolerance.\n",m3);
#endif
                  continue;
                }
                
                // discard this star if matching triangle center is too far from reference triangle center
                // check pointing error along X axis
                d=(ast->x[m1] + ast->x[m2] + ast->x[m3]) / 3;
                if((d2=fabs(d - cenX)) > error){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for excesive X pointing error: %f.\n",m3,d2);
#endif
                  continue;
                }
                // check pointing error along Y axis
                d=(ast->y[m1] + ast->y[m2] + ast->y[m3]) / 3;
                if((d2=fabs(d - cenY)) > error){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for excesive Y pointing error: %f.\n",m3,d2);
#endif
                  continue;
                }
                
                // scale was fixed when found 2nd match candidate, no different scale error is allowed for 3rd matching star 
                // discard if distance between m3 and m1 is not very similar to distance between s3 and s1
                d=sqrt(pow(ast->x[m3]-ast->x[m1],2)+pow(ast->y[m3]-ast->y[m1],2));                
                if(fabs(d - dis13) > DISTMATCHTOL){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for distance to 1st match star error: %f against %f.\n",m3,r,d);
#endif
                  continue;
                }
                // discard if distance between m3 and m2 is not very similar to distance between s3 and s2
                if(fabs(sqrt(pow(ast->x[m3]-ast->x[m2],2)+pow(ast->y[m3]-ast->y[m2],2)) - dis23) > DISTMATCHTOL){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for distance to 2nd match star error.\n",m3);
#endif
                  continue;
                }
                // check angle
                d2=atan((ast->y[m3]-ast->y[m1]) / (ast->x[m3]-ast->x[m1]));
                if((ast->x[m3]-ast->x[m1]) < 0){
                  if(d2 >= 0) d2 -= M_PI;
                  else d2 += M_PI;
                }
                if(fabs(d2 - ang13) > MAXROT){
#if DEBUG>2
                  fprintf(stdout,"Match 3rd star %d discarded for rotation error.\n",m3);
#endif
                  continue;
                }
                // m3 is viable candidate for match to 3rd star
                r=ast->mag[m3];
#if DEBUG>1
                fprintf(stdout,"Match found with 3rd match star: %d, magnitude: %.1f \n",m3,r);
#endif
#if DEBUG>0
			          fprintf(stdout,"Reference triangle: %d, %d, %d\n",s1,s2,s3);
                fprintf(stdout,"Match found with stars: %d, %d, %d\n",m1,m2,m3);
#endif
                triast[2]=m3;
                
                // ************************************************************************************************************
                // At this point a viable candidate triangle is found to match the selected reference triangle
                // will try to validate this candidate by calculating the tranformation between the ref triangle and this candidate
                // We will need to calculate scale, rotation and translation
                
                end_t=clock();
#if DEBUG>1
                fprintf(stdout,"\n");
                fprintf(stdout,"Match triangle found: %d,%d,%d\n",m1,m2,m3); 
                fprintf(stdout,"Magnitudes: %f, %f, %f\n",ast->mag[m1],ast->mag[m2],ast->mag[m3]); 
                fprintf(stdout,"Reference triangle: %d, %d, %d\n",s1,s2,s3);
                fprintf(stdout,"Magnitudes: %f, %f, %f\n",cat->mag[s1],cat->mag[s2],cat->mag[s3]);          
                                          
                d=(float)(end_t - start_t)/CLOCKS_PER_SEC;
                fprintf(stdout,"Seconds elapsed: %f\n",d);          
#endif                        
                start_t=clock();

                // Scale transformation calculation
                // scaling from ref 1-2 to match 1-2 was previously calculated in matchscale, now we average with
                // ref 2-3 to match 2-3 and ref 3-1 to match 3-1
                r=sqrt(pow(cat->x[s3]-cat->x[s2],2)+pow(cat->y[s3]-cat->y[s2],2));
                d=sqrt(pow(ast->x[m3]-ast->x[m2],2)+pow(ast->y[m3]-ast->y[m2],2)) / r;
                matchscale += d;
                r=sqrt(pow(cat->x[s1]-cat->x[s3],2)+pow(cat->y[s1]-cat->y[s3],2));
                d=sqrt(pow(ast->x[m1]-ast->x[m3],2)+pow(ast->y[m1]-ast->y[m3],2)) / r;
                matchscale += d;
                matchscale /= 3;
                
                // Rotation calculation
                // matchrot was obtained from angle between ref stars 1-2 and match stars 1-2
                // now average with angle ref 2-3 to match 2-3 and angle ref 3-1 to match 3-1
                d2=atan((ast->y[m2]-ast->y[m1]) / (ast->x[m2]-ast->x[m1]));
                if((ast->x[m2]-ast->x[m1]) < 0){
                  if(d2 >= 0) d2 -= M_PI;
                  else d2 += M_PI;
                }
                matchrot = d2-ang12;
                d=atan((cat->y[s3]-cat->y[s2]) / (cat->x[s3]-cat->x[s2]));
                if((cat->x[s3]-cat->x[s2]) < 0){
                  if(d >= 0) d -= M_PI;
                  else d += M_PI;
                }
                d2=atan((ast->y[m3]-ast->y[m2]) / (ast->x[m3]-ast->x[m2]));
                if((ast->x[m3]-ast->x[m2]) < 0){
                  if(d2 >= 0) d2 -= M_PI;
                  else d2 += M_PI;
                }
                matchrot += d2-d;
                d=atan((cat->y[s1]-cat->y[s3]) / (cat->x[s1]-cat->x[s3]));
                if((cat->x[s1]-cat->x[s3]) < 0){
                  if(d >= 0) d -= M_PI;
                  else d += M_PI;
                }
                d2=atan((ast->y[m1]-ast->y[m3]) / (ast->x[m1]-ast->x[m3]));
                if((ast->x[m1]-ast->x[m3]) < 0){
                  if(d2 >= 0) d2 -= M_PI;
                  else d2 += M_PI;
                }
                matchrot += d2-d;
                matchrot /= 3;
                
                // traslation calculation
                // we apply averaged matchscale and matchrotation to reference stars and then
                // average reference 1 to match 1 vector, ref2 to match2 vector and ref3 to match3 vector
                r=sqrt(pow(cat->x[s1],2)+pow(cat->y[s1],2));
                d=atan(cat->y[s1] / cat->x[s1]);
                r *= matchscale;
                d += matchrot;
                matchXtras = ast->x[m1] - r*cos(d);
                matchYtras = ast->y[m1] - r*sin(d);
                r=sqrt(pow(cat->x[s2],2)+pow(cat->y[s2],2));
                d=atan(cat->y[s2] / cat->x[s2]);
                r *= matchscale;
                d += matchrot;
                matchXtras += ast->x[m2] - r*cos(d);
                matchYtras += ast->y[m2] - r*sin(d);
                r=sqrt(pow(cat->x[s3],2)+pow(cat->y[s3],2));
                d=atan(cat->y[s3] / cat->x[s3]);
                r *= matchscale;
                d += matchrot;
                matchXtras += ast->x[m3] - r*cos(d);
                matchYtras += ast->y[m3] - r*sin(d);
                matchXtras /= 3;
                matchYtras /= 3;

#if DEBUG>1
                  fprintf(stdout,"Transformation:\nScale: %f\nRotation: %f\nXTraslat: %f\nYTraslat: %f\n",matchscale,matchrot,matchXtras,matchYtras);
#endif
                
                // At this point we've determined the transformation function to convert our reference triangle in the image
                // to a matching triangle in the astrometric catalog
                // Now it's time to check if this tranformation effectively maps our imaged stars catalog into the astrometric catalog
                // For this we apply the tranformation to every bright star in the imaged catalog and accumulate success count when
                // a matching star is found in the astrometric catalog. This count is later used for validation.
                n=0;
                n2=0;
//                d2=(mmin+mave)/2;
//                d2=(mmax+mave)/2;
                d2=mave;
                for(i=0;i<=cat->n;i++){
                  // skip faint stars from the match check
                  if(cat->mag[i] > d2) continue;
                  r=sqrt(pow(cat->x[i],2)+pow(cat->y[i],2));
                  d=atan(cat->y[i] / cat->x[i]);
                  r *= matchscale;
                  d += matchrot;
                  x=r*cos(d) + matchXtras;
                  y=r*sin(d) + matchYtras;
                  n2++;                     // this star counts as a checked star
                  for(j=0;j<=ast->n;j++){
                    // skip also faint stars from map
                    if(ast->mag[j] > (d2 + MAXMAGERR)) continue;  
                    if((fabs(ast->x[j] - x) < DISTMATCHTOL) && (fabs(ast->y[j] - y) < DISTMATCHTOL)){
                      n++;
                      break;
                    }
                  }
                }
                if(n < n2*matchvalratio){
#if DEBUG>0
                  fprintf(stdout,"Match discarded for insufficient match number: %d out of %d checked stars\n",n,n2);
#endif
//                  if(m1==101 && m2==80 && m3==87){
//                    end_t=clock();
//                    fprintf(stdout,"Special match triangle was discarded, halting...\n");
//                    nmatch=3;
//                  }

                  continue;
                }
                fprintf(stdout,"Match count: %d out of %d checked stars\n",n,n2);
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
            // if no match was found for this ref triangle continue looking for ref triangles
            // next try with following 3rd star
            *nselect=2;
#if DEBUG>1
            fprintf(stdout,"No match found for this triangle\n");
#endif
          }
        }
        // if no triangle selected we will try next 2nd star and try again all possible 3rd stars
        if(nmatch<3){
          tricat[2]=0;
          *nselect=1;
#if DEBUG>2
          fprintf(stdout,"Try next 3rd star.\n");
#endif
        }
      }
      // if no triangle selected we will try next 1st star and try again all possible 2nd stars
      if(nmatch<3){
        tricat[1]=0;
        *nselect=0;
#if DEBUG>2
        fprintf(stdout,"Try next 2nd star.\n");
#endif
      }
    }
    // at this point either found 3 candidate stars or already tried all candidate stars without success
    // if loop continues will look for smaller triangles, accept higher magnitude stars and
    // accept lower match ratio for next iteration
    if(nmatch<3){
      size/=1.25;
      mave+=0.25;
    	// relax match constraint
      matchvalratio/=1.25;  
      tricat[0]=0;
#if DEBUG>0
      fprintf(stdout,"Catalog completely parsed with no suitable reference triangle found. Decrasing minimum size, increase max magnitude for next search.\n");
#endif
    }
  }
  // if we could not find suitable triplet erase triange in catalog tricat[] 
  if(nmatch<3){
    tricat[0]=0;
    tricat[1]=0;
    tricat[2]=0;
#if DEBUG>0
    fprintf(stdout,"Did not find any suitable reference triangle.\n");
#endif
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
  
	// Adjust magnitude limit if configured for so
	if (AUTOMAGLIM){
		while (ast.n > cat.n){
			mag -= 0.5;
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
    }
		while (ast.n < cat.n){
			mag += 0.25;
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
    }
    fprintf(stdout,"Astrometric map magnitude limit: %2.2f\n",mag);            
    fprintf(stdout,"Stars in map: %d, extracted imaged stars: %d\n",ast.n,cat.n);
	}

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
#if DEBUG>1
        fprintf(stdout,"\nRef.triangle: ");
        fprintf(stdout,"%d, ",tricat[0]);
#endif
        for(i=1;(i<3 && i<nselect);i++){
          cpgdraw(cat.x[tricat[i]],cat.y[tricat[i]]);
#if DEBUG>1
          fprintf(stdout,"%d, ",tricat[i]);
#endif
        }
        cpgdraw(cat.x[tricat[0]],cat.y[tricat[0]]);
        // astrometric catalog
        if(triast[2]>0){
          cpgmove(ast.x[triast[0]],ast.y[triast[0]]);
#if DEBUG>1
          fprintf(stdout,"\nMatch triangle: ");
          fprintf(stdout,"%d, ",triast[0]);
#endif
          for(i=1;(i<3 && i<nselect);i++){
            cpgdraw(ast.x[triast[i]],ast.y[triast[i]]);
#if DEBUG>1
            fprintf(stdout,"%d, ",triast[i]);
#endif
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
      for (i=0;i<cat.n;i++) {
        cat.select[i]=0;
      }
      for (i=0;i<ast.n;i++) {
        ast.select[i]=0;
      }
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
//      for(i=0;i<3;i++){
//        tricat[i]=0;
//        triast[i]=0;
//        img.a[i]=0;
//        img.b[i]=0;
//      }
//      triast[2]=-1;
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
      //ast=reread_astrometric_catalog(starfile,mag);
      matchvalratio=DEFMATCHVALRATIO;
			tricat[0]=0;
			tricat[1]=0;
			tricat[2]=0;
			triast[0]=0;
			triast[1]=0;
			triast[2]=0;
			
      redraw=1;
      nselect=0;
    }

    // Select pixel catalog
    if (c=='a' && click==0) {
      i=select_nearest(cat,x,y);
      fprintf(stdout,"%d\n",i);      
      cat.select[i]=nselect+1;
      if(nselect<4) 
        tricat[nselect]=i;
      redraw=1;
      click=1;
    }

    // Select catalog
    if (c=='b' && click==1) {
      i=select_nearest(ast,x,y);
      fprintf(stdout,"%d\n",i);            
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
      fprintf(stdout,"Imaged star %d, Mag. %f\n",i,cat.mag[i]);      
    }

    // Select catalog
    if (c=='2') {
      i=select_nearest(ast,x,y);
      fprintf(stdout,"Asto.Catalog star %d, Mag. %f\n",i,ast.mag[i]);            
    }

    // Autoidentify triangles in imaged and astrometric catalogs
    if (c=='i') {
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
      identify_triangles(&cat,&ast,&nselect,tricat,triast);
      // if match was found mark stars into selected fields of catalogs
      if((nselect >=3) && (triast[2] > 0)){
        // erase prior selection first
        for (i=0;i<cat.n;i++) {
          cat.select[i]=0;
        }
        for (i=0;i<ast.n;i++) {
          ast.select[i]=0;
        }
        cat.select[tricat[0]]=1;        
        cat.select[tricat[1]]=2;
        cat.select[tricat[2]]=3;
        ast.select[triast[0]]=1;
        ast.select[triast[1]]=2;
        ast.select[triast[2]]=3;
        nselect=3;
      }
      redraw=1;
    }

    // Reset reference stars selection
    if (c=='I' && nselect>=3) {
      tricat[0]=0;
      triast[0]=0;
      tricat[1]=0;
      triast[1]=0;
      tricat[2]=0;
      triast[2]=-1;
      nselect=0;
      redraw=1;
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
    if (c=='z') {
      width/=1.25;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // Unzoom
    if (c=='x') {
      width*=1.25;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width;
      ymax=y+0.5*width;
      redraw=1;
      continue;
    }

    // More astrometric map stars
    if (c=='+' || c=='=') {
      mag += 0.25;
      fprintf(stdout,"Astrometric map magnitude limit: %2.2f\n",mag);            
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
      fprintf(stdout,"Stars in map: %d, extracted imaged stars: %d\n",ast.n,cat.n);
      redraw=1;
      continue;
    }

    // Less astrometric map stars
    if (c=='-') {
      mag -= 0.25;
      fprintf(stdout,"Astrometric map magnitude limit: %2.2f\n",mag);            
      ast=read_astrometric_catalog(starfile,mag,sx,sy,-q);
      fprintf(stdout,"Stars in map: %d, extracted imaged stars: %d\n",ast.n,cat.n);
      redraw=1;
      continue;
    }

    // Help
    if (c=='h') {
      printf("Calibrates astrometry.\nRequires matching of at least three stars. Use 'i' repeatedly for automatic selection of matching sets or 'a' to manually select a star from the image and 'b' to select the corresponding star from the catalog (tree sets needed). Then use 'f' to check fit, 'R' to start again.\nFinish with 'm' to match stars and 'q' writes calibration output\n");
      printf("q     Quit\n");
      printf("a     Select star on image\n");
      printf("b     Select star from catalog\n");
      printf("i     Autoselect calibration stars from catalog\n");
      printf("I     Reset calibration stars selection\n");
      printf("c     Center image on pixel\n");
      printf("f     Fit calibration\n");
      printf("m     Match stars using current calibration\n");
      printf("z/x   Zoom in/out on cursor\n");
      printf("+/-   Increase/decrease astrometric map stars\n");
      printf("p     Plot sextractor catalog\n");
      printf("r     Reset zoom\n");
      printf("R     Reset fit\n");      
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
