#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qfits.h"
#include "cel.h"
#include "cpgplot.h"

#define LIM 80
#define NMAX 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI

struct fourframe {
  char filename[64];
  int naxis1,naxis2,naxis3,nframes;
  float *zavg,*zstd,*zmax,*znum,*ztrk;
  int *mask;
  double ra0,de0;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  double mjd;
  float *dt,exptime;
  char nfd[32];
  int cospar;
};
struct image {
  int nx,ny;
  float *z;
  float xmin,xmax,ymin,ymax;
};
struct observation {
  int satno,cospar;
  char desig[16],conditions,behavior;
  double mjd,ra,de;
  float terr,perr,tmid;
  char nfd[32],pos[32];
  int epoch,type;
  char iod_line[80];
  float x[3],y[3];
  int state;
};

struct fourframe read_fits(char *filename);

struct image hough_transform(struct image img,int nt,int nr)
{
  int i,j,k,l,ii,jj;
  struct image h;
  float x,y;

  // Limits
  h.xmin=0.0;
  h.xmax=M_PI;
  h.ymin=-sqrt(img.nx*img.nx+img.ny*img.ny);
  h.ymax=sqrt(img.nx*img.nx+img.ny*img.ny);

  // Initialize
  h.nx=nt;
  h.ny=nr;
  h.z=(float *) malloc(h.nx*h.ny*sizeof(float));
  for (i=0;i<h.nx*h.ny;i++)
    h.z[i]=0.0;

  // Loop over image
  for (i=0;i<img.nx;i++) {
    for (j=0;j<img.ny;j++) {
      k=i+img.nx*j;

      if (img.z[k]<1.0)
	continue;

      for (l=0;l<h.nx;l++) {
	x=M_PI*(float) l/(float) (h.nx-1);
	y=(float) i*cos(x)+(float) j*sin(x);
	ii=(int) ((x-h.xmin)/(h.xmax-h.xmin)*h.nx);
	jj=(int) ((y-h.ymin)/(h.ymax-h.ymin)*h.ny);

	if (ii>0 && ii<h.nx && jj>0 && jj<h.ny)
	  h.z[ii+h.nx*jj]+=1.0;
      }
    }
  }

  return h;
}

float hough_peak(struct image h,float *a)
{
  int i,j,k,imax,jmax;
  float hmax;
  float x,y;

  for (i=0;i<h.nx;i++) {
    for (j=0;j<h.ny;j++) {
      k=i+h.nx*j;
      if (k==0 || h.z[k]>hmax) {
	imax=i;
	jmax=j;
	hmax=h.z[k];
      }
    }
  }
  x=h.xmin+(h.xmax-h.xmin)*(float) imax/(float) (h.nx-1);
  y=h.ymin+(h.ymax-h.ymin)*(float) jmax/(float) (h.ny-1);
  a[0]=y/sin(x);
  a[1]=-cos(x)/sin(x);

  return hmax;
}

void update_mask(struct fourframe *ff,float *a,float w,int l)
{
  int i,j,k;
  float phi,sa,ca;
  float x,y,x0,y0;

  // Angle
  phi=-atan2(a[1],1.0);
  sa=sin(phi);
  ca=cos(phi);

  // Loop over pixels
  for (i=0;i<ff->naxis1;i++) {
    for (j=0;j<ff->naxis2;j++) {
      k=i+ff->naxis1*j;
      x0=(float) i;
      y0=(float) j-a[0];
      x=ca*x0-sa*y0;
      y=sa*x0+ca*y0;
      if (fabs(y)<0.5*w && ff->mask[k]==1)
	ff->mask[k]=l;
    }
  }

  return;
}

// Linear least squares fit                                                     
float linear_fit(float x[],float y[],float w[],int n,float a[],float sa[])
{
  int i;
  float sum,sumx,sumy,sumxx,sumxy;
  float d,chi2,covar,r;

  // Compute sums                                                               
  sum=sumx=sumy=sumxx=sumxy=0.;
  for (i=0;i<n;i++) {
    sum+=w[i];
    sumx+=x[i]*w[i];
    sumy+=y[i]*w[i];
    sumxx+=x[i]*x[i]*w[i];
    sumxy+=x[i]*y[i]*w[i];
  }
  d=sum*sumxx-sumx*sumx;

  // Parameters                                                                 
  a[0]=(sumxx*sumy-sumx*sumxy)/d;
  a[1]=(sum*sumxy-sumx*sumy)/d;

  // Uncertainties                                                              
  sa[0]=sqrt(sumxx/d);
  sa[1]=sqrt(sum/d);

  // Chi squared                                                                
  for (i=0,chi2=0.0;i<n;i++)
    chi2+=pow(y[i]-a[0]-a[1]*x[i],2);

  // Covariance                                                                 
  covar= -sumx/d;

  // Correlation coefficient                                                    
  r= -sumx/sqrt(sum*sumxx);

  return chi2;
}

// Get a x and y from a RA and Decl
void forward(double ra0,double de0,double ra,double de,double *x,double *y)
{
  int i;
  char pcode[4]="TAN";
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
  *x*=3600.;
  *y*=3600.;

  return;
}

// Get a RA and Decl from x and y
void reverse(double ra0,double de0,double x,double y,double *ra,double *de)
{
  int i;
  char pcode[4]="TAN";
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
  double mjd,dday;
  float sec;

  sscanf(date,"'%04d-%02d-%02dT%02d:%02d:%f'",&year,&month,&day,&hour,&min,&sec);

  dday=day+hour/24.0+min/1440.0+sec/86400.0;
  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Compute Date from Julian Day
void mjd2date(double mjd,char *date)
{
  double f,jd,dday;
  int z,alpha,a,b,c,d,e;
  double year,month,day,hour,min;
  double sec,x,fsec;

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
  fsec=1000.0*(sec-floor(sec));
  sprintf(date,"%04d%02d%02d%02d%02d%02.0f%03.0f",(int) year,(int) month,(int) day,(int) hour,(int) min,floor(sec),fsec);

  return;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int type)
{
  int i;
  double sec,deg,min,fmin;
  char sign;

  sign=(x<0 ? '-' : '+');
  x=60.*fabs(x);

  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;
  if (type==0)
    fmin=1000.0*(min-floor(min));
  else
    fmin=100.0*(min-floor(min));

  if (type==0)
    sprintf(s,"%02.0f%02.0f%03.0f",deg,floor(min),fmin);
  else
    sprintf(s,"%c%02.0f%02.0f%02.0f",sign,deg,floor(min),fmin);

  return;
}

// Reduce point
void reduce_point(struct observation *obs,struct fourframe img,float tmid,float x,float y)
{
  int iframe,k;
  double ra,de,rx,ry;
  float dx,dy,dt;
  double mjd;
  char nfd[32],sra[15],sde[15];

  // Transform position
  dx=x-img.x0;
  dy=y-img.y0;
  rx=img.a[0]+img.a[1]*dx+img.a[2]*dy;
  ry=img.b[0]+img.b[1]*dx+img.b[2]*dy;
  reverse(img.ra0,img.de0,rx,ry,&ra,&de);

  dec2sex(ra/15.0,sra,0);
  dec2sex(de,sde,1);

  // Get time
  k=(int) x + img.naxis1*(int) y;
  iframe=(int) img.znum[k];
  if (tmid<0.0)
    dt=img.dt[iframe];
  else
    dt=tmid;
  mjd=nfd2mjd(img.nfd)+(double) dt/86400.0;
  mjd2date(mjd,nfd);
  
  // Copy
  strcpy(obs->nfd,nfd);
  sprintf(obs->pos,"%s%s",sra,sde);

  return;
}

void fit(struct observation *obs,struct fourframe ff,int mask)
{
  int i,j,k,l,n,m;
  float *t,*dt,*x,*y,*w;
  float tmin,tmax,tmid;
  float chi2x,chi2y,ax[2],sax[2],ay[2],say[2];
  float dx,dy,dr,rmsx,rmsy;

  // Iterate
  for (m=0;m<10;m++) {
    // Count number of points
    for (i=0,n=0;i<ff.naxis1*ff.naxis2;i++)
      if (ff.mask[i]==mask)
	n++;
    
    // Allocate
    t=(float *) malloc(sizeof(float)*n);
    dt=(float *) malloc(sizeof(float)*n);
    x=(float *) malloc(sizeof(float)*n);
    y=(float *) malloc(sizeof(float)*n);
    w=(float *) malloc(sizeof(float)*n);
    
    // Fill
    for (i=0,l=0;i<ff.naxis1;i++) {
      for (j=0;j<ff.naxis2;j++) {
	k=i+ff.naxis1*j;
	if (ff.mask[k]==mask) {
	  x[l]=(float) i+0.5;
	  y[l]=(float) j+0.5;
	  w[l]=(ff.zmax[k]-ff.zavg[k])/ff.zstd[k];
	  t[l]=ff.dt[(int) ff.znum[k]];
	  l++;
	}
      }
    }
    
    // Find limits in time
    for (i=0;i<n;i++) {
      if (i==0) {
	tmin=t[i];
	tmax=t[i];
      } else {
	if (t[i]<tmin) tmin=t[i];
	if (t[i]>tmax) tmax=t[i];
      }
    }
    tmid=0.5*(tmin+tmax);
    
    // Shift in time
    for (i=0;i<n;i++)
      dt[i]=t[i]-tmid;
    
    // Fit x-pixel position
    chi2x=linear_fit(dt,x,w,n,ax,sax);
    
    // Fit x-pixel position
    chi2y=linear_fit(dt,y,w,n,ay,say);
    
    // Compute rms
    for (i=0,rmsx=0.0,rmsy=0.0;i<n;i++) {
      rmsx+=pow(x[i]-(ax[0]+ax[1]*dt[i]),2);
      rmsy+=pow(y[i]-(ay[0]+ay[1]*dt[i]),2);
    }
    rmsx=sqrt(rmsx/(float) (n-1));
    rmsy=sqrt(rmsy/(float) (n-1));

    // Update mask
    for (i=0;i<ff.naxis1;i++) {
      for (j=0;j<ff.naxis2;j++) {
	k=i+ff.naxis1*j;
	if (ff.mask[k]==mask) {
	  dx=(float) i+0.5-(ax[0]+ax[1]*(ff.dt[(int) ff.znum[k]]-tmid));
	  dy=(float) j+0.5-(ay[0]+ay[1]*(ff.dt[(int) ff.znum[k]]-tmid));
	  if (dx>3.0*rmsx || dy>3.0*rmsy)
	    ff.mask[k]=0;
	}
      }
    }
  }
  printf("%d points, rms (%g, %g), chi2 (%g, %g)\n",n,rmsx,rmsy,chi2x,chi2y);
  //  printf("x: %6.2f +- %.2f  %7.3f +- %.3f  %8.2f; %f pix/s\n",ax[0],sax[0],ax[1],sax[1],chi2x,ax[1]/ff.nframes*ff.exptime);
  //  printf("y: %6.2f +- %.2f  %7.3f +- %.3f  %8.2f; %f pix/s\n",ay[0],say[0],ay[1],say[1],chi2y,ay[1]/ff.nframes*ff.exptime);

  obs->x[0]=ax[0];
  obs->y[0]=ay[0];
  obs->x[1]=ax[0]+ax[1]*(tmin-tmid);
  obs->y[1]=ay[0]+ay[1]*(tmin-tmid);
  obs->x[2]=ax[0]+ax[1]*(tmax-tmid);
  obs->y[2]=ay[0]+ay[1]*(tmax-tmid);
  obs->state=1;

  // Reduce point
  reduce_point(obs,ff,tmid,ax[0],ay[0]);

  // Free
  free(t);
  free(dt);
  free(x);
  free(y);

  return;
}

void format_iod_line(struct observation *obs)
{
  int mt,xt,mp,xp;
  char string[10];
  
  // Time format
  sprintf(string,"%7.1e",obs->terr);
  mt=string[0]-'0';
  xt=atoi(string+4)+8;

  // Position format
  if (obs->type==2) {
    sprintf(string,"%7.1e",obs->perr);
    mp=string[0]-'0';
    xp=atoi(string+4)+8;
  } else {
    printf("Position format not implemented!\n");
  }

  sprintf(obs->iod_line,"%05d %c%c %-6s %04d %c %-17s %d%d %d%d %-14s %d%d %c",
	  obs->satno,
	  obs->desig[0],obs->desig[1],
	  obs->desig+2,
	  obs->cospar,
	  obs->conditions,
	  obs->nfd,
	  mt,xt,
	  obs->type,obs->epoch,
	  obs->pos,
	  mp,xp,
	  obs->behavior);

  return;
}

void write_observation(struct observation obs)
{
  FILE *file;

  file=fopen("detect.dat","a");
  fprintf(file,"%s\n",obs.iod_line);
  fclose(file);

  return;
}

int main(int argc,char *argv[])
{
  int i,j,k,imax,jmax,kmax;
  int flag=0;
  struct fourframe ff;
  struct image img,mask,hough;
  float sigma=5.0,hmax,a[2],hall,fraq;
  struct observation obs;
  char *env;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float zavg,zstd,zmin,zmax;
  char filename[LIM],text[128];

  env=getenv("ST_COSPAR");
  // Default observation
  obs.satno=99999;
  strcpy(obs.desig,"99999U");
  obs.cospar=atoi(env);
  obs.conditions='G';
  strcpy(obs.nfd,"YYYYMMDDHHMMSSsss");
  obs.terr=0.1;
  strcpy(obs.pos,"HHMMmmm+DDMMmm");
  strcpy(obs.iod_line,"");
  obs.perr=0.3;
  obs.epoch=5;
  obs.type=2;
  obs.behavior='S';
  obs.state=0;

  // Read
  if (argc==2)
    ff=read_fits(argv[1]);
  else
    ff=read_fits("test1.fits");

  // Determine limits
  for (i=0,zavg=0.0;i<ff.naxis1*ff.naxis2;i++)
    zavg+=ff.zmax[i];
  zavg/=(float) ff.naxis1*ff.naxis2;
  for (i=0,zstd=0.0;i<ff.naxis1*ff.naxis2;i++)
    zstd+=pow(ff.zmax[i]-zavg,2);
  zstd=sqrt(zstd/(float) (ff.naxis1*ff.naxis2));
  zmin=zavg-2*zstd;
  zmax=zavg+6*zstd;

  // Create image
  img.nx=ff.naxis1;
  img.ny=ff.naxis2;
  img.z=(float *) malloc(img.nx*img.ny*sizeof(float));
  for (i=0;i<ff.naxis1*ff.naxis2;i++) 
    img.z[i]=(ff.zmax[i]-ff.zavg[i])/ff.zstd[i]>sigma;

  // Loop over lines
  for (k=0;;k++) {
    // Generate mask
    for (i=0,hall=0;i<img.nx*img.ny;i++) {
      if ((ff.zmax[i]-ff.zavg[i])/ff.zstd[i]>sigma && ff.mask[i]>=0) {
	img.z[i]=1.0;
	ff.mask[i]=1;
	hall++;
      } else {
	img.z[i]=0.0;
      }
    }

    // Compute Hough transform
    hough=hough_transform(img,1000,1000);

    // Find highest peak
    hmax=hough_peak(hough,a);

    // End
    if (hmax/hall<0.01 || hmax<20)
      break;

    if (flag==0) {
      sprintf(filename,"%s.png/png",argv[1]);
      cpgopen(filename);
      cpgpap(0.,1.0);
      cpgsvp(0.1,0.95,0.1,0.8);
      
      cpgsch(0.8);
      sprintf(text,"UT Date: %.23s  COSPAR ID: %04d",ff.nfd+1,ff.cospar);
      cpgmtxt("T",6.0,0.0,0.0,text);
      sprintf(text,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')",ff.ra0,ff.xrms,ff.de0,ff.yrms);
      cpgmtxt("T",4.8,0.0,0.0,text);
      sprintf(text,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d",ff.naxis1*sqrt(ff.a[1]*ff.a[1]+ff.b[1]*ff.b[1])/3600.0,ff.naxis2*sqrt(ff.a[2]*ff.a[2]+ff.b[2]*ff.b[2])/3600.0,sqrt(ff.a[1]*ff.a[1]+ff.b[1]*ff.b[1]),sqrt(ff.a[2]*ff.a[2]+ff.b[2]*ff.b[2]));
      cpgmtxt("T",3.6,0.0,0.0,text);
      sprintf(text,"Stat: %5.1f+-%.1f (%.1f-%.1f)",zavg,zstd,zmin,zmax);
      cpgmtxt("T",2.4,0.0,0.0,text);
      
      cpgsch(1.0);
      cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      cpgwnad(0.0,(float) ff.naxis1,0.0,(float) ff.naxis2);
      cpgimag(ff.zmax,ff.naxis1,ff.naxis2,1,ff.naxis1,1,ff.naxis2,zmin,zmax,tr);
      //cpgwnad(0.0,(float) img.nx,0.0,(float) img.ny);
      //  cpgimag(img.z,img.nx,img.ny,1,img.nx,1,img.ny,0.0,5.0,tr);
      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);
      flag=1;
    }
      
    printf("Hmax: %.0f/%.0f (%g), %g + %g x\n",hmax,hall,hmax/hall,a[0],a[1]);
  
    // Update mask
    update_mask(&ff,a,20.0,2);

    // Fit
    fit(&obs,ff,2);

    // Format observation
    format_iod_line(&obs);

    // Write observation
    write_observation(obs);

    // Plot observation
    cpgsci(4);
    cpgpt1(obs.x[0],obs.y[0],4);
    cpgmove(obs.x[1],obs.y[1]);
    cpgdraw(obs.x[2],obs.y[2]);
    cpgsci(1);

    // Reset mask
    for (i=0;i<ff.naxis1*ff.naxis2;i++) 
      if (ff.mask[i]>1)
	ff.mask[i]=-1;
  }
  if (flag==1)
    cpgend();

  // Free
  free(img.z);
  free(hough.z);
  free(ff.zavg);
  free(ff.zstd);
  free(ff.zmax);
  free(ff.znum);
  free(ff.dt);
  free(ff.mask);

  return 0;
}

// Read fits fourframe
struct fourframe read_fits(char *filename)
{
  int i,j,k,l,m;
  qfitsloader ql;
  char key[FITS_LINESZ+1];
  char val[FITS_LINESZ+1];
  struct fourframe img;

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
    img.dt[i]=atof(qfits_query_hdr(filename,key));
  }

  // Allocate image memory
  img.zavg=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.zstd=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.zmax=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.znum=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  if (img.naxis3==5) 
    img.ztrk=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.mask=(int *) malloc(sizeof(int)*img.naxis1*img.naxis2);

  // Set mask
  for (i=0;i<img.naxis1*img.naxis2;i++)
    img.mask[i]=0;

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

