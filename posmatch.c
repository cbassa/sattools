#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sgdp4h.h"
#include "satutl.h"
#include <getopt.h>

#define LIM 80

// Cross product
xyz_t cross(xyz_t a,xyz_t b)
{
  xyz_t c;

  c.x=a.y*b.z-a.z*b.y;
  c.y=a.z*b.x-a.x*b.z;
  c.z=a.x*b.y-a.y*b.x;

  return c;
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

// DOY to MJD
double doy2mjd(int year,double doy)
{
  int month,k=2;
  double day;

  if (year%4==0 && year%400!=0)
    k=1;

  month=floor(9.0*(k+doy)/275.0+0.98);
  
  if (doy<32)
    month=1;

  day=doy-floor(275.0*month/9.0)+k*floor((month+9.0)/12.0)+30.0;

  return date2mjd(year,month,day);
}

int main(int argc,char *argv[])
{
  int i,arg=0,imode,n=172800,satno,satname,satno2,imin,flag,status;
  char catalog1[]="catalog.txt",catalog2[]="ison.txt";
  FILE *file1,*file2;
  orbit_t orb1,orb2;
  double mjd0=56908.0,mjd;
  xyz_t *r,*v,r0,v0,d,n1,n2;
  double dr,dv,dn,drmin,dvmin,mjdmin;
  char line0[LIM],line1[LIM],line2[LIM];

  // Reroute stderr
  freopen("/tmp/stderr.txt","w",stderr);

  // Allocate
  r=(xyz_t *) malloc(sizeof(xyz_t)*n);
  v=(xyz_t *) malloc(sizeof(xyz_t)*n);

  // Open file
  file1=fopen(catalog1,"r");
  file2=fopen(catalog2,"r");

  // Loop over classfd catalog
  while (read_twoline(file1,0,&orb1)==0) {
    // Compute positions
    imode=init_sgdp4(&orb1);
    for (i=0;i<n;i++) 
      satpos_xyz(mjd0+2400000.5+(double) i/86400.0,&r[i],&v[i]);

    // Loop over ISON catalog
    flag=0;
    rewind(file2);
    while (read_twoline(file2,0,&orb2)==0) {
      mjd=doy2mjd(orb2.ep_year,orb2.ep_day);

      imode=init_sgdp4(&orb2);
      satpos_xyz(mjd+2400000.5,&r0,&v0);

      // Find nearest offset
      for (i=0;i<n;i++) {
	d.x=r[i].x-r0.x;
	d.y=r[i].y-r0.y;
	d.z=r[i].z-r0.z;
	dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	d.x=v[i].x-v0.x;
	d.y=v[i].y-v0.y;
	d.z=v[i].z-v0.z;
	dv=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);

	if (flag==0 || dr<drmin) {
	  drmin=dr;
	  dvmin=dv;
	  satno=orb2.satno;
	  imin=i;
	  mjdmin=mjd;
	  flag=1;
	}   
      }
    }

    // Compute normals
    rewind(file2);
    read_twoline(file2,satno,&orb2);
    mjd=doy2mjd(orb2.ep_year,orb2.ep_day);
    imode=init_sgdp4(&orb2);
    satpos_xyz(mjd+2400000.5,&r0,&v0);
    n1=cross(r[0],v[0]);
    n2=cross(r0,v0);
    d.x=n1.x-n2.x;
    d.y=n1.y-n2.y;
    d.z=n1.z-n2.z;
    dn=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);

    // Find object name
    rewind(file2);
    while (fgetline(file2,line0,LIM)>0) {
      fgetline(file2,line1,LIM);
      fgetline(file2,line2,LIM);
      sscanf(line1,"1 %d",&satno2);
      if (satno==satno2) {
	sscanf(line0,"SO %d",&satname);
      }
    }

    //    printf("%05d %05d %8.2f km %9.5f km/s %6.0lf s %8.0f km\n",orb1.satno,satno,drmin,dvmin,(mjd0+(double) imin/86400.0-mjd)*86400,dn);
    printf("%6d %05d %s | %5d %8.2f km %9.5f km/s %6.0lf s %8.0f km\n",satname,orb1.satno,orb1.desig,satno,drmin,dvmin,(mjd0+(double) imin/86400.0-mjd)*86400,dn);
  }
  fclose(file1);
  fclose(file2);

  // Free
  free(r);
  free(v);

  return 0;
}
