#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <wcslib/cel.h>
#include <cpgplot.h>
#include "qfits.h"

#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define LIM 256

struct image {
  char filename[64];
  int naxis1,naxis2,naxis3,nframes;
  float *zavg,*zstd,*zmax,*znum,*zd,*zsig;
  int *mask;
  char nfd[32];
  double ra0,de0;
  float x0,y0;
  float a[3],b[3],xrms,yrms;
  double mjd;
  float *dt,exptime;
  int cospar;
};
struct selection {
  int state,fit;
  float x0,y0,x1,y1;
  float w,zmin,zmax;
  float a,ca,sa,r;
  float tmid,tmin,tmax,ax[2],sax[2],ay[2],say[2],chi2x,chi2y;
};
struct observation {
  int satno,cospar;
  char desig[16],conditions,behavior;
  double mjd,ra,de;
  float terr,perr,tmid,tmin,tmax;
  char nfd[32],pos[32];
  int epoch,type;
  char iod_line[80];
  float x[3],y[3];
  float ax[2],ay[2];
  int state;
};
struct track {
  float x0,y0,x1,y1,texp;
  int satno;
} trk;
int iobject=0;
  
struct image read_fits(char *filename);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
void reverse(double ra0,double de0,double x,double y,double *ra,double *de);
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void mjd2date(double mjd,char *date);
void dec2sex(double x,char *s,int type);
float linear_fit(float x[],float y[],int n,float a[],float sa[]);
int fgetline(FILE *file,char *s,int lim);
double gmst(double mjd);
double modulo(double x,double y);

// Find peak
float find_peak(float *z,int kx,int ky,int xmin,int xmax,int ymin,int ymax,float s,int mx,int my,float *x0,float *y0)
{
  int i,j,k,l,i0,j0,k0,i1,j1,k1,nx,ny;
  int imid,jmid,imax,jmax,flag;
  float *d,*g,*w,*h;
  float sg,sgg,sgn,den,s1,s2;
  float hmax,havg,hstd;
  float x,y;

  printf("%d %d %d %d -> %d %d\n",xmin,xmax,ymin,ymax,kx,ky);
  // Select image
  if (xmin<0.0) xmin=0.0;
  if (ymin<0.0) ymin=0.0;
  if (xmax>=kx) xmax=kx;
  if (ymax>=ky) ymax=ky;
  nx=(int) (xmax-xmin);
  ny=(int) (ymax-ymin);
  d=(float *) malloc(sizeof(float)*nx*ny);
  for (i=xmin,i0=0;i<xmax;i++,i0++) {
    for (j=ymin,j0=0;j<ymax;j++,j0++) {
      k=i+kx*j;
      k0=i0+nx*j0;
      d[k0]=z[k];
    }
  }
  printf("%d %d %d %d -> %d %d\n",xmin,xmax,ymin,ymax,nx,ny);
  // Make kernel
  g=(float *) malloc(sizeof(float)*mx*my);
  for (i=0,imid=mx/2,jmid=my/2,sg=0.0,sgg=0.0;i<mx;i++) {
    x=((float) (i-imid))/s;
    for (j=0;j<my;j++) {
      y=((float) (j-jmid))/s;
      k=i+mx*j;
      g[k]=exp(-0.5*(x*x+y*y));
      sg+=g[k];
      sgg+=g[k]*g[k];
    }
  }
  sgn=sg/(float) (mx*my);
  den=sgg-sg*sg/(float) (mx*my);

  // Compute weights
  w=(float *) malloc(sizeof(float)*mx*my);
  for (i=0;i<mx;i++) {
    for (j=0;j<my;j++) {
      k=i+mx*j;
      w[k]=(g[k]-sgn)/den;
    }
  }

  // Compute H array
  h=(float *) malloc(sizeof(float)*nx*ny);
  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      k=i+nx*j;
      h[k]=0.0;
      for (i0=0;i0<mx;i0++) {
	for (j0=0;j0<my;j0++) {
	  k0=i0+mx*j0;
	  i1=i+i0-imid;
	  j1=j+j0-jmid;
	  // Keep in range
	  if (i1<0 || i1>=nx || j1<0 || j1>=ny)
	    continue;
	  k1=i1+nx*j1;
	  h[k]+=w[k0]*d[k1];
	}
      }
      h[k]/=(float) (mx*my);
    }
  }
  
  // Locate maximum
  for (i=mx,flag=0;i<nx-mx;i++) {
    for (j=my;j<ny-my;j++) {
      k=i+nx*j;
      if (flag==0) {
	imax=i;
	jmax=j;
	hmax=h[k];
	flag=1;
      }
      if (h[k]>hmax) {
	imax=i;
	jmax=j;
	hmax=h[k];
      }
    }
  }

  // Compute mean value
  for (i=mx,s1=0.0,s2=0.0;i<nx-mx;i++) {
    for (j=my;j<ny-my;j++) {
      k=i+nx*j;
      s1+=h[k];
      s2+=1.0;
    }
  }
  havg=s1/s2;

  // Standard deviation
  for (i=mx,s1=0.0,s2=0.0;i<nx-mx;i++) {
    for (j=my;j<ny-my;j++) {
      k=i+nx*j;
      s1+=(h[k]-havg)*(h[k]-havg);
      s2+=1.0;
    }
  }
  hstd=sqrt(s1/s2);

  // Free
  free(g);
  free(w);
  free(h);
  free(d);

  *x0=imax+xmin+0.5;
  *y0=jmax+ymin+0.5;

  return (hmax-havg)/hstd;
}

// MJD to DOY
double mjd2doy(double mjd,int *yr)
{
  int year,month,k=2;
  int day;
  double doy;
  char nfd[32];
  
  mjd2date(mjd,nfd);

  sscanf(nfd,"%04d",&year);
  sscanf(nfd+4,"%02d",&month);
  sscanf(nfd+6,"%02d",&day);

  if (year%4==0 && year%400!=0)
    k=1;

  doy=floor(275.0*month/9.0)-k*floor((month+9.0)/12.0)+day-30;

  *yr=year;

  return doy;
}

// Reduce point
void reduce_point(struct observation *obs,struct image img,float tmid,float x,float y)
{
  int i,iframe,k;
  double ra,de,rx,ry;
  float dx,dy,dt;
  double mjd,mjd1,dra;
  char nfd[32],sra[15],sde[15];
  float ax[2],ay[2];
  
  // Transform position
  dx=x-img.x0;
  dy=y-img.y0;
  rx=img.a[0]+img.a[1]*dx+img.a[2]*dy;
  ry=img.b[0]+img.b[1]*dx+img.b[2]*dy;
  reverse(img.ra0,img.de0,rx,ry,&ra,&de);

  // Transform direction
  for (i=0;i<2;i++) {
    ax[i]=obs->ax[i];
    ay[i]=obs->ay[i];
  }
  obs->ax[1]=(img.a[1]*ax[1]+img.a[2]*ay[1])/3600.0;
  obs->ay[1]=(img.b[1]*ax[1]+img.b[2]*ay[1])/3600.0;
  
  // Get time
  k=(int) x + img.naxis1*(int) y;
  iframe=(int) img.znum[k];
  if (tmid<0.0)
    dt=img.dt[iframe];
  else
    dt=tmid;
  mjd=nfd2mjd(img.nfd)+(double) dt/86400.0;
  mjd2date(mjd,nfd);
  obs->mjd=mjd;

  // Correct for motion
  mjd1=img.mjd+0.5*(double) img.exptime/86400.0;
  dra=gmst(mjd)-gmst(mjd1);
  ra+=dra;

  // Get RA/Dec
  dec2sex(ra/15.0,sra,0);
  dec2sex(de,sde,1);
  obs->ra=ra;
  obs->de=de;
  
  // Copy
  strcpy(obs->nfd,nfd);
  sprintf(obs->pos,"%s%s",sra,sde);

  return;
}

void compute_cuts(float *z,int *mask,int n,float *zmin,float *zmax,float lcut,float hcut)
{
  int i,m;
  double s1,s2;
  double avg,std;

  for (i=0,s1=0.0,s2=0.0,m=0;i<n;i++) {
    if (mask[i]==1) {
      s1+=z[i];
      s2+=z[i]*z[i];
      m++;
    }
  }
  avg=s1/(float) m;
  std=sqrt(s2/(float) m-avg*avg);
  
  *zmin=(float) (avg-lcut*std);
  *zmax=(float) (avg+hcut*std);
  
  return;
}

void plot_selection(struct selection s)
{
  int i;
  float x,y,dx,dy;

  cpgsci(7);
  if (s.state==1) {
    cpgpt1(s.x0,s.y0,17);
  } 
  //  if (s.state>1) {
  //    cpgmove(s.x0,s.y0);
  //    cpgdraw(s.x1,s.y1);
  //  }

  if (s.state==2) {
    for (i=0;i<5;i++) {
      if (i==0 || i==4) {
	dx=-s.w;
	dy=-s.w;
      } else if (i==1) {
	dx=-s.w;
	dy=s.w;
      } else if (i==2) {
	dx=s.w;
	dy=s.w;
      } else if (i==3) {
	dx=s.w;
	dy=-s.w;
      }
      dx=0.0;
      if (i<2 || i==4) {
	x=s.ca*dx-s.sa*dy+s.x0;
	y=s.sa*dx+s.ca*dy+s.y0;
      } else {
	x=s.ca*dx-s.sa*dy+s.x1;
	y=s.sa*dx+s.ca*dy+s.y1;
      }
      
      if (i==0)
	cpgmove(x,y);
      else
	cpgdraw(x,y);

    }
  }
  cpgsci(1);



  return;
}

void apply_mask(struct image *img,struct selection s)
{
  int i,j,k;
  float x,y,dx,dy;

  for (i=0;i<img->naxis1;i++) {
    for (j=0;j<img->naxis2;j++) {
      k=i+img->naxis1*j;
      if (img->mask[k]==0)
	continue;
      dx=(float) i-s.x0;
      dy=(float) j-s.y0;
      x=s.ca*dx+s.sa*dy;
      y=-s.sa*dx+s.ca*dy;
      if (x>=0.0 && x<=s.r && y>-s.w && y<s.w && img->zmax[k]>s.zmin) {
	img->mask[k]=1;
      } else {
	img->mask[k]=0;
      }
    }
  }
  
  return;
}

void apply_mask_sigma(struct image *img,struct selection s)
{
  int i,j,k;
  float x,y,dx,dy;

  for (i=0;i<img->naxis1;i++) {
    for (j=0;j<img->naxis2;j++) {
      k=i+img->naxis1*j;
      if (img->mask[k]==0)
	continue;
      dx=(float) i-s.x0;
      dy=(float) j-s.y0;
      x=s.ca*dx+s.sa*dy;
      y=-s.sa*dx+s.ca*dy;
      if (x>=0.0 && x<=s.r && y>-s.w && y<s.w && img->zsig[k]>s.zmin) {
	img->mask[k]=1;
      } else {
	img->mask[k]=0;
      }
    }
  }
  
  return;
}

void mask_pixel(struct image *img,float x,float y)
{
  int i,j,k,kmin,i0,j0,flag;
  float r,rmin;

  i0=(int) x;
  j0=(int) y;

  // Find nearest pixel
  for (i=0,flag=0;i<img->naxis1;i++) {
    for (j=0;j<img->naxis2;j++) {
      k=i+img->naxis1*j;
      r=sqrt(pow(i-i0,2)+pow(j-j0,2));
      if (img->mask[k]==0)
	continue;
      if (flag==0 || r<rmin) {
	rmin=r;
	kmin=k;
	flag=1;
      }
    }
  }

  // Mask pixel
  img->mask[kmin]=0;

  return;
}



void fit(struct observation *obs,struct image img)
{
  int i,j,k,l,n;
  float *t,*dt,*x,*y;
  float tmin,tmax,tmid;
  float chi2x,chi2y,ax[2],sax[2],ay[2],say[2];

  // Count number of points
  for (i=0,n=0;i<img.naxis1*img.naxis2;i++)
    if (img.mask[i]==1)
      n++;

  // Allocate
  t=(float *) malloc(sizeof(float)*n);
  dt=(float *) malloc(sizeof(float)*n);
  x=(float *) malloc(sizeof(float)*n);
  y=(float *) malloc(sizeof(float)*n);

  // Fill
  for (i=0,l=0;i<img.naxis1;i++) {
    for (j=0;j<img.naxis2;j++) {
      k=i+img.naxis1*j;
      if (img.mask[k]==1) {
	x[l]=(float) i+0.5;
	y[l]=(float) j+0.5;
	t[l]=img.dt[(int) img.znum[k]];
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
  printf("Using points between %.3f and %.3f\n",tmin,tmax);
    
  // Shift in time
  for (i=0;i<n;i++)
    dt[i]=t[i]-tmid;

  // Fit x-pixel position
  chi2x=linear_fit(dt,x,n,ax,sax);

  // Fit x-pixel position
  chi2y=linear_fit(dt,y,n,ay,say);

  printf("x: %6.2f +- %.2f  %7.3f +- %.3f  %8.2f; %f pix/s\n",ax[0],sax[0],ax[1],sax[1],chi2x,ax[1]/img.nframes*img.exptime);
  printf("y: %6.2f +- %.2f  %7.3f +- %.3f  %8.2f; %f pix/s\n",ay[0],say[0],ay[1],say[1],chi2y,ay[1]/img.nframes*img.exptime);

  obs->x[0]=ax[0];
  obs->y[0]=ay[0];
  obs->x[1]=ax[0]+ax[1]*(tmin-tmid);
  obs->y[1]=ay[0]+ay[1]*(tmin-tmid);
  obs->x[2]=ax[0]+ax[1]*(tmax-tmid);
  obs->y[2]=ay[0]+ay[1]*(tmax-tmid);
  obs->state=1;
  obs->tmin=tmin;
  obs->tmax=tmax;
  for (i=0;i<2;i++) {
    obs->ax[i]=ax[i];
    obs->ay[i]=ay[i];
  }
  
  // Reduce point
  reduce_point(obs,img,tmid,ax[0],ay[0]);

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

void find_designation(int satno0,char *desig0)
{
  FILE *file;
  int satno;
  char desig[16];
  char *env,filename[128];

  // Environment variables
  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/desig.txt",env);

  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"Designation file not found!\n");
    exit(0);
  }
  while (!feof(file)) {
    fscanf(file,"%d %s",&satno,desig);
    if (satno==satno0) {
      strcpy(desig0,desig);
      break;
    }
  }
  fclose(file);

  return;
}

void write_observation(struct observation obs)
{
  FILE *file;
  float w,pa,dt;
  
  file=fopen("observations.txt","a");
  fprintf(file,"%s\n",obs.iod_line);
  fclose(file);

  printf("Observation written\n");

  dt=obs.tmax-obs.tmin;
  w=sqrt(obs.ax[1]*obs.ax[1]+obs.ay[1]*obs.ay[1])/dt;
  pa=atan2(obs.ay[1],obs.ax[1])*R2D;
  
  return;
}

void track(char *fileroot,struct observation obs,struct image *img,float frac)
{
  FILE *file;
  char line[LIM],filename[LIM];
  int flag=0,satno;
  float x0,y0,x1,y1,texp;
  int i,j,k,l,k0;
  int di,dj;
  float *z;
  int *wt;
  float dxdn,dydn,dx,dy;

  sprintf(filename,"%s.id",fileroot);

  // Open ID file
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"ID file %s not found\n",filename);
    return;
  }
  while (fgetline(file,line,LIM)>0) {
    sscanf(line,"%s %f %f %f %f %f %d",filename,&x0,&y0,&x1,&y1,&texp,&satno);
    trk.x0=x0;
    trk.y0=y0;
    trk.x1=x1;
    trk.y1=y1;
    trk.satno=satno;
    trk.texp=texp;
    if (satno==obs.satno)
      break;
  }
  fclose(file);

  if (satno!=obs.satno) {
    fprintf(stderr,"Object %d not found\n",obs.satno);
    return;
  }
  dxdn=(x1-x0)/(float) img->nframes;
  dydn=(y1-y0)/(float) img->nframes;

  // Allocate
  z=(float *) malloc(sizeof(float)*img->naxis1*img->naxis2);
  wt=(int *) malloc(sizeof(int)*img->naxis1*img->naxis2);

  // Set to zero
  for (i=0;i<img->naxis1*img->naxis2;i++) {
    z[i]=0.0;
    wt[i]=0;
  }

  // Loop over frames
  for (l=0;l<img->nframes;l++) {
    // Offset
    dx=dxdn*(l-frac*img->nframes);
    dy=dydn*(l-frac*img->nframes);
    
    // Integer offset
    di=(int) floor(dx+0.5);
    dj=(int) floor(dy+0.5);

    // Set
    for (i=0;i<img->naxis1;i++) {
      for (j=0;j<img->naxis2;j++) {
	k=i+img->naxis1*j;
	k0=i+di+img->naxis1*(j+dj);
	if (i+di>0 && i+di<img->naxis1 && j+dj>0 && j+dj<img->naxis2) {
	  wt[k]+=1;
	  if (img->znum[k0]==l)
	    z[k]+=img->zmax[k0];
	  //	  else
	  //	    z[k]+=img->zavg[k0];
	}
      }
    }
  }

  // Scale
  for (i=0;i<img->naxis1*img->naxis2;i++) {
    if (wt[i]>0)
      img->zd[i]=z[i]/(float) wt[i];
    else
      img->zd[i]=z[i];
  }
  img->naxis3=5;

  free(z);
  free(wt);

  return;
}

int autotrack(char *fileroot,struct observation obs,struct image *img,int cflag)
{
  FILE *file;
  char line[LIM],filename[LIM];
  int flag=0,satno,satno0=0,i=0,n;
  float x0,y0,x1,y1,texp;
  int status=0;

  sprintf(filename,"%s.id",fileroot);

  // Open ID file
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"ID file %s not found\n",filename);
    return -1;
  }
  while (fgetline(file,line,LIM)>0) {
    if (cflag==1 && strstr(line,"classfd")==NULL)
      continue;
    sscanf(line,"%s %f %f %f %f %f %d",filename,&trk.x0,&trk.y0,&trk.x1,&trk.y1,&trk.texp,&trk.satno);
    if (i==iobject) {
      status=1;
      break;
    }

    i++;
  }
  fclose(file);

  iobject++;

  return status;
}

int main(int argc,char *argv[])
{
  int i,j,k,l;
  int iconditions=0,ibehavior=0;
  struct image img;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float x,y,frac=0.5;
  char c;
  float xmin,xmax,ymin,ymax,zmin,zmax,*z;
  float width;
  int redraw=1,layer=2,status;
  float lcut=4,hcut=6;
  struct selection s;
  struct observation obs;
  char conditions[]="EGFPBT",behavior[]="EFIRSX";
  char text[128];
  double doy,mjd;
  int year;
  char *env;
  float sigma,xa,ya;

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

  // Set track
  trk.satno=0;

  // Read image
  img=read_fits(argv[1]);

  // Allocate
  z=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);

  // Get fake designation
  mjd=nfd2mjd(img.nfd);
  doy=mjd2doy(mjd,&year);
  sprintf(obs.desig,"%02d%03.0lfA",year-2000,doy+500);

  cpgopen("/xs");
  cpgpap(0.,1.0);
  //cpgpap(7,0.75);
  cpgask(0);
  cpgsch(0.8);

  // Default limits
  xmin=0.0;
  xmax=(float) img.naxis1;
  ymin=0.0;
  ymax=(float) img.naxis2;
  width=img.naxis1;

  // Default selection
  s.state=0;
  s.w=10;
  s.zmin=5.0;
  s.zmax=20.0;
  s.fit=0;

  // Set cospas
  obs.cospar=img.cospar;

  for (;;) {
    if (redraw==1) {
      cpgeras();
      
      cpgsvp(0.1,0.95,0.1,0.95);
      cpgwnad(xmin,xmax,ymin,ymax);
      cpglab("x (pix)","y (pix)"," ");
      cpgsfs(2);
      cpgctab (heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);

      sprintf(text,"UT Date: %.23s  COSPAR ID: %04d",img.nfd+1,img.cospar);
      cpgmtxt("T",6.0,0.0,0.0,text);
      sprintf(text,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')",img.ra0,img.xrms,img.de0,img.yrms);
      cpgmtxt("T",4.8,0.0,0.0,text);
      sprintf(text,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d",img.naxis1*sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1])/3600.0,img.naxis2*sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2])/3600.0,sqrt(img.a[1]*img.a[1]+img.b[1]*img.b[1]),sqrt(img.a[2]*img.a[2]+img.b[2]*img.b[2]));
      cpgmtxt("T",3.6,0.0,0.0,text);
  
      // Apply mask
      for (i=0;i<img.naxis1*img.naxis2;i++) {
	if (layer==2) z[i]=img.zmax[i]*img.mask[i];
	if (layer==3) z[i]=img.znum[i]*img.mask[i];
	if (layer==4) z[i]=img.zd[i]*img.mask[i];
	if (layer==5) z[i]=img.zsig[i]*img.mask[i];
      }
      
      if (layer==0) compute_cuts(img.zavg,img.mask,img.naxis1*img.naxis2,&zmin,&zmax,lcut,hcut);
      if (layer==1) compute_cuts(img.zstd,img.mask,img.naxis1*img.naxis2,&zmin,&zmax,lcut,hcut);
      if (layer==2) compute_cuts(img.zmax,img.mask,img.naxis1*img.naxis2,&zmin,&zmax,lcut,hcut);
      if (layer==4) compute_cuts(img.zd,img.mask,img.naxis1*img.naxis2,&zmin,&zmax,lcut,hcut);
      if (layer==5) {
	zmin=s.zmin;
	zmax=s.zmax;
      }

      if (layer==0) cpgimag(img.zavg,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      if (layer==1) cpgimag(img.zstd,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      if (layer==2) cpgimag(z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      if (layer==3) cpgimag(z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,0.0,(float) img.nframes,tr);
      if (layer==4) cpgimag(z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);
      if (layer==5) cpgimag(z,img.naxis1,img.naxis2,1,img.naxis1,1,img.naxis2,zmin,zmax,tr);

      cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);

      // Plot fit
      if (obs.state==1) {
	cpgsci(4);
	cpgpt1(obs.x[0],obs.y[0],4);
	cpgmove(obs.x[1],obs.y[1]);
	cpgdraw(obs.x[2],obs.y[2]);
	cpgsci(1);
      } else if (obs.state==2) {
	cpgsci(4);
	cpgpt1(obs.x[0],obs.y[0],4);
	cpgsci(1);
      }

      // Plot selection
      if (s.state!=0) 
	plot_selection(s);

      // Plot track
      if (trk.satno!=0) {
	cpgsci(4);
	cpgmove(trk.x0,trk.y0);
	cpgdraw(trk.x1,trk.y1);
	cpgpt1(frac*(trk.x1-trk.x0)+trk.x0,frac*(trk.y1-trk.y0)+trk.y0,4);
	cpgsci(1);
      }

      format_iod_line(&obs);
      cpgmtxt("T",1.0,0.5,0.5,obs.iod_line);
      redraw=0;
    }

    

    // Get cursor
    cpgcurs(&x,&y,&c);

    // Quit
    if (c=='q')
      break;

    // End
    if (c=='a') {
      status=autotrack(argv[1],obs,&img,0);
      if (status==1) {
	obs.satno=trk.satno;
	find_designation(obs.satno,obs.desig);
	track(argv[1],obs,&img,frac);
	x=frac*(trk.x1-trk.x0)+trk.x0;
	y=frac*(trk.y1-trk.y0)+trk.y0;
	width=100;
	xmin=x-0.5*width;
	xmax=x+0.5*width;
	ymin=y-0.5*width*img.naxis2/img.naxis1;
	ymax=y+0.5*width*img.naxis2/img.naxis1;
	sigma=find_peak(img.zd,img.naxis1,img.naxis2,(int) xmin,(int) xmax,(int) ymin,(int) ymax,2.0,11,11,&x,&y);
	printf("%f %f %f\n",x,y,sigma);
	if (sigma>5.0) {
	  reduce_point(&obs,img,frac*img.exptime,x,y);
	  obs.x[0]=x;
	  obs.y[0]=y;
	  obs.state=2;
	}

	redraw=1;
	layer=4;
      }
      continue;
    }

    // Find peak
    if (c=='p') {
      sigma=find_peak(img.zd,img.naxis1,img.naxis2,(int) xmin,(int) xmax,(int) ymin,(int) ymax,2.0,11,11,&x,&y);
      printf("%f %f %f\n",x,y,sigma);
      reduce_point(&obs,img,frac*img.exptime,x,y);
      obs.x[0]=x;
      obs.y[0]=y;
      obs.state=2;
      redraw=1;
      continue;
    }
    
    // Track
    if (c=='t') {
      printf("Provide satellite ID: ");
      scanf("%d",&obs.satno);
      find_designation(obs.satno,obs.desig);
      track(argv[1],obs,&img,frac);
      layer=4;
      redraw=1;
    }
	
    if (c=='\t') {
      status=autotrack(argv[1],obs,&img,1);
      if (status==1) {
	obs.satno=trk.satno;
	find_designation(obs.satno,obs.desig);
	track(argv[1],obs,&img,frac);
	redraw=1;
	layer=4;
      }
    }

    // Write obs
    if (c=='w') {
      write_observation(obs);
      continue;
    }

    // Reduce
    if (c=='m') {
      reduce_point(&obs,img,-1.0,x,y);
      obs.x[0]=x;
      obs.y[0]=y;
      obs.state=2;
      redraw=1;
      continue;
    }

    // Change fraction
    if (c=='e') {
      if (frac>0.49 && frac<0.51)
	frac=1.0;
      else if (frac>0.51)
	frac=0.0;
      else if (frac<0.49)
	frac=0.5;
      printf("Fraction: %.1f\n",frac);
      iobject=0;
    }

    // Change fraction
    if (c=='E') {
      frac+=0.1;
      if (frac>1.0)
	frac=0.0;
      printf("Fraction: %.1f\n",frac);
      iobject=0;
    }

    // Reduce
    if (c=='M' || c=='D') {
      reduce_point(&obs,img,frac*img.exptime,x,y);
      obs.x[0]=x;
      obs.y[0]=y;
      obs.state=2;
      redraw=1;
      continue;
    }

    // Get designation
    if (c=='d') {
      printf("Provide satellite number: ");
      scanf("%d",&obs.satno);
      find_designation(obs.satno,obs.desig);
      redraw=1;
      continue;
    }

    // Toggle condition
    if (c=='C') {
      iconditions++;
      if (iconditions>strlen(conditions)-1)
	iconditions=0;
      obs.conditions=conditions[iconditions];
      redraw=1;
      continue;
    }
    // Toggle behavior
    if (c=='B') {
      ibehavior++;
      if (ibehavior>strlen(behavior)-1)
	ibehavior=0;
      obs.behavior=behavior[ibehavior];
      redraw=1;
      continue;
    }

    // Reread
    if (c=='R') {
      img=read_fits(argv[1]);
      redraw=1;
      continue;
    }

    // Start
    if (c=='s' && s.state==0) {
      s.x0=x;
      s.y0=y;
      s.state=1;
      redraw=1;
      continue;
    }

    // Fit
    if (c=='F') {
      fit(&obs,img);
      redraw=1;
      continue;
    }

    // End
    if (c=='f' && s.state==1) {
      s.x1=x;
      s.y1=y;
      s.a=atan2(s.y1-s.y0,s.x1-s.x0);
      s.ca=cos(s.a);
      s.sa=sin(s.a);
      s.r=sqrt(pow(s.x0-s.x1,2)+pow(s.y0-s.y1,2));
      s.state=2;
      apply_mask_sigma(&img,s);
      //s.zmin=zmin;
      redraw=1;
      continue;
    }

    // Mask pixel
    if (c=='X' && s.state!=0) {
      mask_pixel(&img,x,y);
      apply_mask(&img,s);
      redraw=1;
      continue;
    }

    // Change level
    if (c=='+' || c=='=') {
      s.zmin*=1.01;
      printf("%.4f\n",s.zmin);
      apply_mask_sigma(&img,s);
      redraw=1;
      continue;
    }
    if (c=='-') {
      s.zmin/=1.01;
      printf("%.4f\n",s.zmin);
      apply_mask_sigma(&img,s);
      redraw=1;
      continue;
    }

    // Mean
    if (isdigit(c)) {
      layer=c-'0'-1;
      redraw=1;
      continue;
    }

    // Adjust cuts
    if (c=='v') {
      lcut*=2;
      hcut*=2;
      redraw=1;
      continue;
    }
    if (c=='b') {
      lcut/=2;
      hcut/=2;
      if (lcut<0.5) lcut=0.5;
      if (hcut<0.75) hcut=0.75;
      redraw=1;
      continue;
    }

    // Center
    if (c=='c') {
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width*img.naxis2/img.naxis1;
      ymax=y+0.5*width*img.naxis2/img.naxis1;
      redraw=1;
      continue;
    }

    // Zoom
    if (c=='z') {
      width/=2;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width*img.naxis2/img.naxis1;
      ymax=y+0.5*width*img.naxis2/img.naxis1;
      redraw=1;
      continue;
    }

    // Unzoom
    if (c=='x') {
      width*=2;
      xmin=x-0.5*width;
      xmax=x+0.5*width;
      ymin=y-0.5*width*img.naxis2/img.naxis1;
      ymax=y+0.5*width*img.naxis2/img.naxis1;
      redraw=1;
      continue;
    }

    // Reset
    if (c=='r') {
      xmin=0.0;
      xmax=(float) img.naxis1;
      ymin=0.0;
      ymax=(float) img.naxis2;
      width=img.naxis1;
      lcut=4.0;
      hcut=6.0;
      s.state=0;
      s.fit=0;
      obs.state=0;
      redraw=1;
      continue;
    }

		if (c=='h'){
      printf("Reduce Satellite tracks. ");
      printf("q   Quit\n");
			printf("TAB track and stack objects automatically (only classfd sats)\n");
			printf("w		write IOD observation to observations.txt\n");
			printf("M/D measure stack and track position (also middle mouse button)\n");
			printf("e   change fraction (0.0 to start, 0.5 for medium, 1.0 for the end)\n"); 
			printf("d   provide NORAD satellite number\n");
			printf("C   toggle IOD observing conditions G-Good F-Fair P-Poor B-Bad T-Terrible E-Excellent\n");
			printf("B   toggle behavior of sat.: F-Flash I-Irregular R-Regular S-Steady X-Uspecified E-Extremely weak\n");
			printf("s   select start of satellite track\n");
			printf("f   select end of satellite track\n");
			printf("F   fit satellite track\n");
			printf("+/= Increase level of masking pixels\n");
			printf("-   Lower level for masking pixels (This does not seem to work)\n");
			printf("1   go to the mean pixel value FITS layer\n");
			printf("2   go to the FITS standard deviation layer\n");
			printf("3   go to the maximum pixel value FITS layer\n");
			printf("4   go to the frame number of the maximum pixel value FITS layer\n");
			printf("5   go to the stack and track layer (only after 't' or 'TAB')\n");
			printf("v   lower dynamic range\n");
			printf("b   increase dynamic range\n");
			printf("c   center on cursor\n");
			printf("z   zoom in at cursor\n");
			printf("x   zoom out at cursor\n");
			printf("R/r reset to start\n");
			printf("m   measure position of pixel\n");
			printf("t   track & stack, give NORAD satellite number\n");
			printf("E   change fraction in tenths\n");
			printf("X   pixel mask (also with mouse scroll wheel)\n");
		}
  }

  cpgend();

  free(img.zavg);
  free(img.zstd);
  free(img.zmax);
  free(img.znum);
  free(img.zd);
  free(img.zsig);

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
  img.nframes=atoi(qfits_query_hdr(filename,"NFRAMES"));

  // MJD
  img.mjd=(double) atof(qfits_query_hdr(filename,"MJD-OBS"));
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

  // Timestamps
  img.dt=(float *) malloc(sizeof(float)*img.nframes);
  for (i=0;i<img.nframes;i++) {
    sprintf(key,"DT%04d",i);
    //strcpy(val,qfits_query_hdr(filename,key));
    //    sscanf(val+1,"%f",&img.dt[i]);
    img.dt[i]=atof(qfits_query_hdr(filename,key));
  }

  // Allocate image memory
  img.zavg=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.zstd=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.zmax=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.znum=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.zd=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  img.mask=(int *) malloc(sizeof(int)*img.naxis1*img.naxis2);
  img.zsig=(float *) malloc(sizeof(float)*img.naxis1*img.naxis2);
  for (i=0;i<img.naxis1*img.naxis2;i++) 
    img.zd[i]=0.0;

  // Set parameters
  ql.xtnum=0;
  ql.ptype=PTYPE_FLOAT;
  ql.filename=filename;

  // Loop over planes
  for (k=0;k<img.naxis3;k++) {
    ql.pnum=k;;

    // Initialize load
    if (qfitsloader_init(&ql) != 0) 
      printf("Error initializing data loading\n");

    // Test load
    if (qfits_loadpix(&ql) != 0) 
      printf("Error loading actual data\n");

    // Fill z array
    for (i=0,l=0;i<img.naxis1;i++) {
      for (j=0;j<img.naxis2;j++) {
	if (k==0) img.zavg[l]=ql.fbuf[l];
	if (k==1) img.zstd[l]=ql.fbuf[l];
	if (k==2) img.zmax[l]=ql.fbuf[l];
	if (k==3) img.znum[l]=ql.fbuf[l];
	if (k==4) img.zd[l]=ql.fbuf[l];
	l++;
      }
    }
  }

  // Compute scaled
  for (i=0;i<img.naxis1*img.naxis2;i++) 
    img.mask[i]=1;

  // Compute sigma
  for (i=0;i<img.naxis1*img.naxis2;i++) 
    img.zsig[i]=(img.zmax[i]-img.zavg[i])/img.zstd[i];

  return img;
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
  fsec=floor(1000.0*(sec-floor(sec)));
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
    fmin=floor(1000.0*(min-floor(min)));
  else
    fmin=floor(100.0*(min-floor(min)));

  if (type==0)
    sprintf(s,"%02.0f%02.0f%03.0f",deg,floor(min),fmin);
  else
    sprintf(s,"%c%02.0f%02.0f%02.0f",sign,deg,floor(min),fmin);

  return;
}

// Linear least squares fit                                                     
float linear_fit(float x[],float y[],int n,float a[],float sa[])
{
  int i;
  float sum,sumx,sumy,sumxx,sumxy;
  float w,d,chi2,covar,r;

  // Compute sums                                                               
  sum=sumx=sumy=sumxx=sumxy=0.;
  for (i=0;i<n;i++) {
    w=1.0;
    sum+=w;
    sumx+=x[i]*w;
    sumy+=y[i]*w;
    sumxx+=x[i]*x[i]*w;
    sumxy+=x[i]*y[i]*w;
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

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\t')
    c=' ';
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}
