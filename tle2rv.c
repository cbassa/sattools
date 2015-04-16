#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sgdp4h.h"
#include "satutl.h"
#include <getopt.h>


#define LIM 128
#define XKE 0.07436680 // Guassian Gravitational Constant
#define XKMPER 6378.135
#define AE 1.0
#define XMNPDA 1440.0
#define R2D 180.0/M_PI
#define D2R M_PI/180.0

extern double SGDP4_jd0;

// Nutation series
static const struct {
  int nm1,nm,nf,nd,nn;
  double s,st;
  double c,ct;
} nut[]={
  {  0,  0,  0,  0,  1, -171996.0, -174.2,  92025.0,    8.9 },
  {  0,  0,  0,  0,  2,    2062.0,    0.2,   -895.0,    0.5 },
  { -2,  0,  2,  0,  1,      46.0,    0.0,    -24.0,    0.0 },
  {  2,  0, -2,  0,  0,      11.0,    0.0,      0.0,    0.0 },
  { -2,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 },
  {  1, -1,  0, -1,  0,      -3.0,    0.0,      0.0,    0.0 },
  {  0, -2,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 },
  {  2,  0, -2,  0,  1,       1.0,    0.0,      0.0,    0.0 },
  {  0,  0,  2, -2,  2,  -13187.0,   -1.6,   5736.0,   -3.1 },
  {  0,  1,  0,  0,  0,    1426.0,   -3.4,     54.0,   -0.1 },
  {  0,  1,  2, -2,  2,    -517.0,    1.2,    224.0,   -0.6 },
  {  0, -1,  2, -2,  2,     217.0,   -0.5,    -95.0,    0.3 },
  {  0,  0,  2, -2,  1,     129.0,    0.1,    -70.0,    0.0 },
  {  2,  0,  0, -2,  0,      48.0,    0.0,      1.0,    0.0 },
  {  0,  0,  2, -2,  0,     -22.0,    0.0,      0.0,    0.0 },
  {  0,  2,  0,  0,  0,      17.0,   -0.1,      0.0,    0.0 },
  {  0,  1,  0,  0,  1,     -15.0,    0.0,      9.0,    0.0 },
  {  0,  2,  2, -2,  2,     -16.0,    0.1,      7.0,    0.0 },
  {  0, -1,  0,  0,  1,     -12.0,    0.0,      6.0,    0.0 },
  { -2,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 },
  {  0, -1,  2, -2,  1,      -5.0,    0.0,      3.0,    0.0 },
  {  2,  0,  0, -2,  1,       4.0,    0.0,     -2.0,    0.0 },
  {  0,  1,  2, -2,  1,       4.0,    0.0,     -2.0,    0.0 },
  {  1,  0,  0, -1,  0,      -4.0,    0.0,      0.0,    0.0 },
  {  2,  1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 },
  {  0,  0, -2,  2,  1,       1.0,    0.0,      0.0,    0.0 },
  {  0,  1, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  1,  0,  0,  2,       1.0,    0.0,      0.0,    0.0 },
  { -1,  0,  0,  1,  1,       1.0,    0.0,      0.0,    0.0 },
  {  0,  1,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  0,  2,  0,  2,   -2274.0,   -0.2,    977.0,   -0.5 },
  {  1,  0,  0,  0,  0,     712.0,    0.1,     -7.0,    0.0 },
  {  0,  0,  2,  0,  1,    -386.0,   -0.4,    200.0,    0.0 },
  {  1,  0,  2,  0,  2,    -301.0,    0.0,    129.0,   -0.1 },
  {  1,  0,  0, -2,  0,    -158.0,    0.0,     -1.0,    0.0 },
  { -1,  0,  2,  0,  2,     123.0,    0.0,    -53.0,    0.0 },
  {  0,  0,  0,  2,  0,      63.0,    0.0,     -2.0,    0.0 },
  {  1,  0,  0,  0,  1,      63.0,    0.1,    -33.0,    0.0 },
  { -1,  0,  0,  0,  1,     -58.0,   -0.1,     32.0,    0.0 },
  { -1,  0,  2,  2,  2,     -59.0,    0.0,     26.0,    0.0 },
  {  1,  0,  2,  0,  1,     -51.0,    0.0,     27.0,    0.0 },
  {  0,  0,  2,  2,  2,     -38.0,    0.0,     16.0,    0.0 },
  {  2,  0,  0,  0,  0,      29.0,    0.0,     -1.0,    0.0 },
  {  1,  0,  2, -2,  2,      29.0,    0.0,    -12.0,    0.0 },
  {  2,  0,  2,  0,  2,     -31.0,    0.0,     13.0,    0.0 },
  {  0,  0,  2,  0,  0,      26.0,    0.0,     -1.0,    0.0 },
  { -1,  0,  2,  0,  1,      21.0,    0.0,    -10.0,    0.0 },
  { -1,  0,  0,  2,  1,      16.0,    0.0,     -8.0,    0.0 },
  {  1,  0,  0, -2,  1,     -13.0,    0.0,      7.0,    0.0 },
  { -1,  0,  2,  2,  1,     -10.0,    0.0,      5.0,    0.0 },
  {  1,  1,  0, -2,  0,      -7.0,    0.0,      0.0,    0.0 },
  {  0,  1,  2,  0,  2,       7.0,    0.0,     -3.0,    0.0 },
  {  0, -1,  2,  0,  2,      -7.0,    0.0,      3.0,    0.0 },
  {  1,  0,  2,  2,  2,      -8.0,    0.0,      3.0,    0.0 },
  {  1,  0,  0,  2,  0,       6.0,    0.0,      0.0,    0.0 },
  {  2,  0,  2, -2,  2,       6.0,    0.0,     -3.0,    0.0 },
  {  0,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 },
  {  0,  0,  2,  2,  1,      -7.0,    0.0,      3.0,    0.0 },
  {  1,  0,  2, -2,  1,       6.0,    0.0,     -3.0,    0.0 },
  {  0,  0,  0, -2,  1,      -5.0,    0.0,      3.0,    0.0 },
  {  1, -1,  0,  0,  0,       5.0,    0.0,      0.0,    0.0 },
  {  2,  0,  2,  0,  1,      -5.0,    0.0,      3.0,    0.0 },
  {  0,  1,  0, -2,  0,      -4.0,    0.0,      0.0,    0.0 },
  {  1,  0, -2,  0,  0,       4.0,    0.0,      0.0,    0.0 },
  {  0,  0,  0,  1,  0,      -4.0,    0.0,      0.0,    0.0 },
  {  1,  1,  0,  0,  0,      -3.0,    0.0,      0.0,    0.0 },
  {  1,  0,  2,  0,  0,       3.0,    0.0,      0.0,    0.0 },
  {  1, -1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 },
  { -1, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 },
  { -2,  0,  0,  0,  1,      -2.0,    0.0,      1.0,    0.0 },
  {  3,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 },
  {  0, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 },
  {  1,  1,  2,  0,  2,       2.0,    0.0,     -1.0,    0.0 },
  { -1,  0,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 },
  {  2,  0,  0,  0,  1,       2.0,    0.0,     -1.0,    0.0 },
  {  1,  0,  0,  0,  2,      -2.0,    0.0,      1.0,    0.0 },
  {  3,  0,  0,  0,  0,       2.0,    0.0,      0.0,    0.0 },
  {  0,  0,  2,  1,  2,       2.0,    0.0,     -1.0,    0.0 },
  { -1,  0,  0,  0,  2,       1.0,    0.0,     -1.0,    0.0 },
  {  1,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 },
  { -2,  0,  2,  2,  2,       1.0,    0.0,     -1.0,    0.0 },
  { -1,  0,  2,  4,  2,      -2.0,    0.0,      1.0,    0.0 },
  {  2,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  1,  1,  2, -2,  2,       1.0,    0.0,     -1.0,    0.0 },
  {  1,  0,  2,  2,  1,      -1.0,    0.0,      1.0,    0.0 },
  { -2,  0,  2,  4,  2,      -1.0,    0.0,      1.0,    0.0 },
  { -1,  0,  4,  0,  2,       1.0,    0.0,      0.0,    0.0 },
  {  1, -1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 },
  {  2,  0,  2, -2,  1,       1.0,    0.0,     -1.0,    0.0 },
  {  2,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0 },
  {  1,  0,  0,  2,  1,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  0,  4, -2,  2,       1.0,    0.0,      0.0,    0.0 },
  {  3,  0,  2, -2,  2,       1.0,    0.0,      0.0,    0.0 },
  {  1,  0,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  1,  2,  0,  1,       1.0,    0.0,      0.0,    0.0 },
  { -1, -1,  0,  2,  1,       1.0,    0.0,      0.0,    0.0 },
  {  0,  0, -2,  0,  1,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  0,  2, -1,  2,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  1,  0,  2,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  1,  0, -2, -2,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  0, -1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0 },
  {  1,  1,  0, -2,  1,      -1.0,    0.0,      0.0,    0.0 },
  {  1,  0, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 },
  {  2,  0,  0,  2,  0,       1.0,    0.0,      0.0,    0.0 },
  {  0,  0,  2,  4,  2,      -1.0,    0.0,      0.0,    0.0 },
  {  0,  1,  0,  1,  0,       1.0,    0.0,      0.0,    0.0 }
};

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Present nfd
void nfd_now(char *s)
{
  time_t rawtime;
  struct tm *ptm;

  // Get UTC time
  time(&rawtime);
  ptm=gmtime(&rawtime);
    
  sprintf(s,"%04d-%02d-%02dT%02d:%02d:%02d",ptm->tm_year+1900,ptm->tm_mon+1,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
  
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
  int year,month,day,hour,min,sec;
  double mjd,dday;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%02d",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;

  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Nutation
void nutation(double mjd,double *dpsi,double *deps,double *eps)
{
  int i;
  double t,t2,t3;
  double d,m,m1,f,n;
  double arg;

  // Julian centuries
  t=(mjd-51544.5)/36525.0;
  t2=t*t;
  t3=t2*t;

  // Angles
  d=modulo(297.85036+445267.111480*t-0.0019142*t2+t3/189474.0,360.0)*D2R;
  m=modulo(357.52772+35999.050340*t-0.0001603*t2-t3/300000.0,360.0)*D2R;
  m1=modulo(134.96298+477198.867398*t+0.0086972*t2+t3/56250.0,360.0)*D2R;
  f=modulo(93.27191+483202.017538*t-0.0036825*t2+t3/327270.0,360.0)*D2R;
  n=modulo(125.04452-1934.136261*t+0.0020708*t2+t3/450000.0,360.0)*D2R;
  
  // Compute sums
  *dpsi=0.0;
  *deps=0.0;
  for (i=0;i<106;i++) {
    arg=nut[i].nd*d+nut[i].nm*m+nut[i].nm1*m1+nut[i].nf*f+nut[i].nn*n;
    *dpsi+=(nut[i].s+nut[i].st*t)*sin(arg);
    *deps+=(nut[i].c+nut[i].ct*t)*cos(arg);
  }
  *dpsi*=0.0001/3600*D2R;
  *deps*=0.0001/3600*D2R;
  *eps=-46.8150*t-0.00059*t2+0.001813*t3;
  *eps=(23.4392911+ *eps/3600.0)*D2R;

  return;
}

// Precession
void precess(double mjd0,double mjd,double *zeta,double *z,double *theta)
{
  double t0,t;

  // Time in centuries
  t0=(mjd0-51544.5)/36525.0;
  t=(mjd-mjd0)/36525.0;

  // Precession angles
  *zeta=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *zeta+=(0.30188-0.000344*t0)*t*t+0.017998*t*t*t;
  *zeta*=D2R/3600.0;
  *z=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *z+=(1.09468+0.000066*t0)*t*t+0.018203*t*t*t;
  *z*=D2R/3600.0;
  *theta=(2004.3109-0.85330*t0-0.000217*t0*t0)*t;
  *theta+=-(0.42665+0.000217*t0)*t*t-0.041833*t*t*t;
  *theta*=D2R/3600.0;
  
  return;
}


// Set identity matrix
void identity_matrix(double a[3][3])
{
  int i,j;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      if (i==j)
	a[i][j]=1.0;
      else
	a[i][j]=0.0;
    }
  }

  return;
}

// Rotate around x-axis
void rotate_x(double phi, double a[3][3])
{
  double s,c,a10,a11,a12,a20,a21,a22;

  s=sin(phi);
  c=cos(phi);

  a10=c*a[1][0]+s*a[2][0];
  a11=c*a[1][1]+s*a[2][1];
  a12=c*a[1][2]+s*a[2][2];
  a20=-s*a[1][0]+c*a[2][0];
  a21=-s*a[1][1]+c*a[2][1];
  a22=-s*a[1][2]+c*a[2][2];
  
  a[1][0]=a10;
  a[1][1]=a11;
  a[1][2]=a12;
  a[2][0]=a20;
  a[2][1]=a21;
  a[2][2]=a22;
  
  return;
}

// Rotate around y-axis
void rotate_y(double phi, double a[3][3])
{
  double s,c,a00,a01,a02,a20,a21,a22;

  s=sin(phi);
  c=cos(phi);

  a00=c*a[0][0]-s*a[2][0];
  a01=c*a[0][1]-s*a[2][1];
  a02=c*a[0][2]-s*a[2][2];
  a20=s*a[0][0]+c*a[2][0];
  a21=s*a[0][1]+c*a[2][1];
  a22=s*a[0][2]+c*a[2][2];
  
  a[0][0]=a00;
  a[0][1]=a01;
  a[0][2]=a02;
  a[2][0]=a20;
  a[2][1]=a21;
  a[2][2]=a22;
  
  return;
}

// Rotate around z-axis
void rotate_z(double phi, double a[3][3])
{
  double s,c,a00,a01,a02,a10,a11,a12;

  s=sin(phi);
  c=cos(phi);

  a00=c*a[0][0]+s*a[1][0];
  a01=c*a[0][1]+s*a[1][1];
  a02=c*a[0][2]+s*a[1][2];
  a10=-s*a[0][0]+c*a[1][0];
  a11=-s*a[0][1]+c*a[1][1];
  a12=-s*a[0][2]+c*a[1][2];

  a[0][0]=a00;
  a[0][1]=a01;
  a[0][2]=a02;
  a[1][0]=a10;
  a[1][1]=a11;
  a[1][2]=a12;

  return;
}

// Matrix multiply
void matrix_multiply(double a[3][3],double b[3][3],double c[3][3])
{
  int i,j,k;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      c[i][j]=0.0;
      for (k=0;k<3;k++) 
	c[i][j]+=a[i][k]*b[k][j];
    }
  }

  return;
}

// Vector multiply
void vector_multiply_in_place(double a[3][3],double b[3])
{
  int i,j,k;
  double c[3];

  for (i=0;i<3;i++) {
    c[i]=0.0;
    for (j=0;j<3;j++) 
      c[i]+=a[i][j]*b[j];
  }
  for (i=0;i<3;i++)
    b[i]=c[i];

  return;
}

// Vector multiply
void vector_multiply(double a[3][3],double b[3],double c[3])
{
  int i,j,k;

  for (i=0;i<3;i++) {
    c[i]=0.0;
    for (j=0;j<3;j++) 
      c[i]+=a[i][j]*b[j];
  }

  return;
}

// Transpose
void matrix_transpose(double a[3][3],double b[3][3])
{
  int i,j;

  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      b[i][j]=a[j][i];

  return;
}

// Dot product
double dot_product(double a[3],double b[3])
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}  

// Magnitude
double magnitude(double a[3])
{
  return sqrt(dot_product(a,a));
}

// Cross product
void cross_product(double a[3],double b[3],double c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];

  return;
}

// ICRS to TEME conversion
void icrs_to_teme(double mjd,double a[3][3])
{
  int i,j;
  double dpsi,deps,eps,z,theta,zeta,h;
  double p[3][3],n[3][3],q[3][3],b[3][3];

  // Precession
  precess(51544.5,mjd,&zeta,&z,&theta);
  identity_matrix(p);
  rotate_z(-zeta,p);
  rotate_y(theta,p);
  rotate_z(-z,p);

  // Nutation
  nutation(mjd,&dpsi,&deps,&eps);
  identity_matrix(n);
  rotate_x(eps,n);
  rotate_z(-dpsi,n);
  rotate_x(-eps-deps,n);

  // Equation of equinoxes
  identity_matrix(q);
  rotate_z(dpsi*cos(eps+deps),q);

  // Multiply matrices (left to right)
  matrix_multiply(q,n,b);
  matrix_multiply(b,p,a);

  return;
}

void usage(void)
{
  printf("tle2rv c:i:t:m:efh\n\n");
  printf("-c   Catalog to load [classfd.tle]\n");
  printf("-i   Object to convert [all]\n");
  printf("-t   Epoch (YYYY-mm-ddThh:mm:ss)\n");
  printf("-m   Epoch (MJD)\n");
  printf("-e   Use TLE epoch\n");
  printf("-j   J2000 output\n");
  printf("-g   GMAT output (J2000)\n");
  printf("-h   This help\n");

  return;
}

int main(int argc,char *argv[])
{
  int imode,arg,satno=0,useepoch=0,format=0,gmat=0;
  FILE *file;
  char *env;
  char tlefile[LIM],nfd[32];
  double mjd;
  xyz_t r,v;
  orbit_t orb;
  double rr[3],vv[3],e[3][3],et[3][3];

  // Get environment variable
  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/classfd.tle",env);

  // Set date
  nfd_now(nfd);
  mjd=nfd2mjd(nfd);

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:t:m:hejg"))!=-1) {
    switch (arg) {

    case 't':
      strcpy(nfd,optarg);
      mjd=nfd2mjd(nfd);
      break;
    
    case 'm':
      mjd=atof(optarg);
      break;

    case 'e':
      useepoch=1;
      break;

    case 'j':
      format=1;
      break;

    case 'g':
      gmat=1;
      break;

    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 'i':
      satno=atoi(optarg);
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

  // Open file
  file=fopen(tlefile,"r");
  while (read_twoline(file,satno,&orb)==0) {
    // Propagate
    imode=init_sgdp4(&orb);

    // Use epoch instead of user supplied date
    if (useepoch==1)
      mjd=SGDP4_jd0-2400000.5;

    // Compute position and velocity
    imode=satpos_xyz(mjd+2400000.5,&r,&v);

    // Output
    if (format==0 && gmat==0)
      printf("%05d %14.8lf %f %f %f %f %f %f TEME\n",orb.satno,mjd,r.x,r.y,r.z,v.x,v.y,v.z);

    // To vectors
    rr[0]=r.x;
    rr[1]=r.y;
    rr[2]=r.z;
    vv[0]=v.x;
    vv[1]=v.y;
    vv[2]=v.z;

    // Matrices
    icrs_to_teme(mjd,e);
    matrix_transpose(e,et);

    // Transform
    vector_multiply_in_place(et,rr);
    vector_multiply_in_place(et,vv);
    
    // Output J2000
    if (format==1 && gmat==0)
      printf("%05d %14.8lf %f %f %f %f %f %f J2000\n",orb.satno,mjd,rr[0],rr[1],rr[2],vv[0],vv[1],vv[2]);

    // GMAT output
    if (gmat==1) {
      printf("UTCModJulian = %14.8lf\n",mjd-29999.5);
      printf("CoordinateSystem = EarthMJ2000Eq\n");
      printf("X = %lf\n",rr[0]);
      printf("Y = %lf\n",rr[1]);
      printf("Z = %lf\n",rr[2]);
      printf("VX = %lf\n",vv[0]);
      printf("VY = %lf\n",vv[1]);
      printf("VZ = %lf\n",vv[2]);
    }
  }
  fclose(file);

  return 0;
}
