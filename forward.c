#include <stdio.h>
#include <string.h>
#include <math.h>
#include <wcslib/cel.h>

// Get a x and y from a RA and Decl
void forward(double ra0,double de0,double ra,double de,double *x,double *y)
{
  int i,status;
  double phi,theta;
  struct celprm cel;

  // Initialize Reference Angles
  celini(&cel);
  cel.ref[0]=ra0;
  cel.ref[1]=de0;
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;
  strcpy(cel.prj.code,"STG");

  if (celset(&cel)) {
    printf("Error in Projection (celset)\n");
    return;
  }
  cels2x(&cel,1,0,1,1,&ra,&de,&phi,&theta,x,y,&status);

  *x *=3600.;
  *y *=3600.;

  return;
}
