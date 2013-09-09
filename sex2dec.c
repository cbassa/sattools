// Compute decimal from sexagesimal
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// Usage: sex2dec [options] hh:mm:ss.ss
int main(int argc,char *argv[])
{
  double x;
  int sign=1;
  float deg,min,sec;
  char t[20],c;

  if (argc==1) {
    printf("Usage: sex2dec -r <hh:mm:ss.sss>\n    -or-    <ddd:mm:ss.sss>\nCompute sexagesimal from decimal input x.\n");
    printf("Options: -r Converts hours into degrees\n");
    
    return -1;
  }

  // Get Sexagesimal string
  strcpy(t,argv[--argc]);
  if (t[0]=='-') sign=-1;

  deg=fabs(atof(strtok(t," :,")));
  min=fabs(atof(strtok(NULL," :,")));
  sec=fabs(atof(strtok(NULL," :,")));
  
  x=(double) deg+(double) min/60.+(double) sec/3600.;

  // Get Options
  while (--argc > 0 && (*++argv)[0] == '-') {
    while (c = *++argv[0]) {
      if (c == 'r')
	x *= 15.;
    }
  }
  printf("%lf\n",sign*x);

  return 0;
}
