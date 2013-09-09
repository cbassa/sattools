#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cel.h"

#define LIM 80

void dec2sex(double,char *,int,int);
double sex2dec(char *);
void reverse(double,double,double,double,double *,double *);
void forward(double,double,double,double,double *,double *);

int main(int argc,char *argv[])
{
  int i;
  double ra1,de1,ra2,de2;
  double rx,ry;
  char sra[15],sde[15];

  if (argc==1) {
    printf("Usage: %s <ra1> <de1> <ra2> <de2>\n",argv[0]);
    printf("       Computes angular offset\n");
    printf("Usage: %s -d <ra1> <de1> <dra> <dde>\n",argv[0]);
    printf("       Applies angular offset\n");
    printf("Usage: %s -x <ra1> <de1> <ra2> <de2>\n",argv[0]);
    printf("       Computes x-offset only\n");
    printf("Usage: %s -y <ra1> <de1> <ra2> <de2>\n",argv[0]);
    printf("       Computes y-offset only\n");
    return -1;
  }

  if (strcmp(argv[1],"-d")==0) {
    if (strchr(argv[2],':')!=NULL)
      ra1=15.*sex2dec(argv[2]);
    else
      ra1=atof(argv[2]);
    if (strchr(argv[3],':')!=NULL)
      de1=sex2dec(argv[3]);
    else
      de1=atof(argv[3]);

    rx=(double) atof(argv[4]);
    ry=(double) atof(argv[5]);

    reverse(ra1,de1,rx,ry,&ra2,&de2);

    dec2sex(ra2/15.0,sra,0,7);
    dec2sex(de2,sde,0,6);

    printf("%s %s\n",sra,sde);
  } else if (strcmp(argv[1],"-x")==0) {
    if (strchr(argv[2],':')!=NULL)
      ra1=15.*sex2dec(argv[2]);
    else
      ra1=atof(argv[2]);
    if (strchr(argv[3],':')!=NULL)
      de1=sex2dec(argv[3]);
    else
      de1=atof(argv[3]);
    if (strchr(argv[4],':')!=NULL)
      ra2=15.*sex2dec(argv[4]);
    else
      ra2=atof(argv[4]);
    if (strchr(argv[5],':')!=NULL)
      de2=sex2dec(argv[5]);
    else
      de2=atof(argv[5]);
    
    forward(ra1,de1,ra2,de2,&rx,&ry);
    
    printf("%8.3f\n",rx);
  } else if (strcmp(argv[1],"-y")==0) {
    if (strchr(argv[2],':')!=NULL)
      ra1=15.*sex2dec(argv[2]);
    else
      ra1=atof(argv[2]);
    if (strchr(argv[3],':')!=NULL)
      de1=sex2dec(argv[3]);
    else
      de1=atof(argv[3]);
    if (strchr(argv[4],':')!=NULL)
      ra2=15.*sex2dec(argv[4]);
    else
      ra2=atof(argv[4]);
    if (strchr(argv[5],':')!=NULL)
      de2=sex2dec(argv[5]);
    else
      de2=atof(argv[5]);
    
    forward(ra1,de1,ra2,de2,&rx,&ry);
    
    printf("%8.3f\n",ry);
  } else {
    if (strchr(argv[1],':')!=NULL)
      ra1=15.*sex2dec(argv[1]);
    else
      ra1=atof(argv[1]);
    if (strchr(argv[2],':')!=NULL)
      de1=sex2dec(argv[2]);
    else
      de1=atof(argv[2]);
    if (strchr(argv[3],':')!=NULL)
      ra2=15.*sex2dec(argv[3]);
    else
      ra2=atof(argv[3]);
    if (strchr(argv[4],':')!=NULL)
      de2=sex2dec(argv[4]);
    else
      de2=atof(argv[4]);
    
    forward(ra1,de1,ra2,de2,&rx,&ry);
    
    printf("dRA:%8.3f dDE:%8.3f d:%8.3f\n",rx,ry,sqrt(rx*rx+ry*ry));
  }

  return 0;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int f,int len)
{
  int i;
  double sec,deg,min;
  char sign;
  char *form[]={":: ",",, ","hms","   "};

  sign=(x<0 ? '-' : ' ');
  x=3600.*fabs(x);

  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;

  if (len==7) sprintf(s,"%c%02i%c%02i%c%07.4f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==6) sprintf(s,"%c%02i%c%02i%c%06.3f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==5) sprintf(s,"%c%02i%c%02i%c%05.2f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==4) sprintf(s,"%c%02i%c%02i%c%04.1f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==2) sprintf(s,"%c%02i%c%02i%c%02i%c",sign,(int) deg,form[f][0],(int) min,form[f][1],(int) floor(sec),form[f][2]);

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
  *x *=3600.;
  *y *=3600.;

  return;
}
