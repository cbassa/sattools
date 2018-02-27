#include <stdio.h>
#include <string.h>
#include <math.h>
#include <wcslib/cel.h>

// Get a RA and Decl from x and y
void reverse(double ra0,double de0,double x,double y,double *ra,double *de)
{
  int i,status;
  double phi,theta;
  struct celprm cel;

  x/=3600.;
  y/=3600.;

  // Initialize Reference Angles
  celini(&cel);
  cel.ref[0]=ra0;
  cel.ref[1]=de0;
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;
  strcpy(cel.prj.code,"TAN");
  
  if (celset(&cel)) {
    printf("Error in Projection (celset)\n");
    return;
  }
  celx2s(&cel,1,0,1,1,&x,&y,&phi,&theta,ra,de,&status);

  return;
}
