#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "qfits.h"

int main(int argc, char * argv[])
{
  int i;
  char keyword[FITS_LINESZ+1];
  char *value;
  
  // Usage
  if (argc<3) {
    printf("Usage: %s <filename> [ext] <key1> <key2> etc.\n", argv[0]);
    return 1 ;
  }

  // Check this is indeed a FITS file 
  if (is_fits_file(argv[1])!=1) {
    printf("%s is not a FITS file\n", argv[1]);
    return -1 ;
  }
  
  // Extension header?
  if (atoi(argv[2])==0) {
    for (i=2;i<argc;i++) 
      printf("%s ",qfits_query_hdr(argv[1], argv[i]));
  } else {
    for (i=3;i<argc;i++)
      printf("%s ",qfits_query_ext(argv[1], argv[i],atoi(argv[2])));
  }
  printf("\n");


  return 0 ;
}
