/* Compute sexagesimal from decimal */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc,char *argv[])
{
  int h=0;
  double sec,deg,min;
  char c,sign,d=' ';
  double x;
  char s[14];

  if (argc==1) {
    printf("Usage: dec2sex [options] <x>\n\nCompute sexagesimal from decimal input x.\n");
    printf("Options: -r gives hours instead of degrees\n");
    printf("         -d<a> uses character a as delimiter\n");
    printf("         -h uses hms as delimiters\n");
    printf("         -s always prints sign\n\n");
    
    return -1;
  }
  
  // Get Decimal Value
  x=(double) atof(argv[--argc]);
  sign=(x<0 ? '-' : ' ');

  // Get Options
  while (--argc > 0 && (*++argv)[0] == '-') {
    while (c = *++argv[0]) {
      if (c == 'r')
	x /= 15.;
      if (c == 's')
	sign=(x<0 ? '-' : '+');
      if (c == 'd')
	if (strlen(argv[0])!=1)
	  d=(*argv)[1];
      if (c == 'h')
	h++;
    }
  }

  x=3600.*fabs(x);
  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  deg=x;

  if (h==0)
    printf("%c%02i%c%02i%c%06.3f\n",sign,(int) deg,d,(int) min,d,(float) sec);
  else
    printf("%c%02ih%02im%06.3fs\n",sign,(int) deg,(int) min,(float) sec);

  return 0;
}
