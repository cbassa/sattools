#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <ctype.h>
#include "sgdp4h.h"
#include "satutl.h"

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;
 
  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  //  if (c == '\n')
  //    s[i++] = c;
  s[i] = '\0';
  return i;
}

void usage(void)
{
  // while ((arg=getopt(argc,argv,"c:i:I:D:N:h"))!=-1) {
  printf("Usage: mvtle -i OLD_NORAD_ID -I NEW_NORAD_ID -d COSPAR_ID [-N SATELLITE_NAME] [-h]\n\n");
  printf("Arguments:\n");
  printf("-i OLD_NORAD_ID     Old Satellite Catalog Number (NORAD ID)\n");
  printf("-i NEW_NORAD_ID     New Satellite Catalog Number (NORAD ID)\n");
  printf("-d COSPAR_ID        International Designator (COSPAR ID)\n");
  printf("-N SATELLITE_NAME   (New) Satellite Name / TLE Line 0, default: ''\n");
  printf("-h                  This help\n");

  return;
}

int main(int argc,char *argv[])
{
  int arg=0;
  char *tlefile;
  int satno=0,satno_new=0;
  char line0[70]="",line1[70],line2[70],desig[10]="";
  orbit_t orb;
  FILE *file;

  // Decode options
  while ((arg=getopt(argc,argv,"c:i:I:D:N:h"))!=-1) {
    switch (arg) {
      
    case 'c':
      tlefile=optarg;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'I':
      satno_new=atoi(optarg);
      break;

    case 'D':
      strcpy(desig,optarg);
      break;

    case 'N':
      strcpy(line0,optarg);
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
  file=fopen(tlefile,"rb");
  if (file==NULL) 
    fatal_error("File open failed for reading \"%s\"",tlefile);
  
  // Loop over file
  while (read_twoline(file,satno,&orb)==0) {
    // Adjust changes
    if (satno_new!=0) 
      orb.satno=satno_new;
    if (strcmp(desig,"")!=0) 
      strcpy(orb.desig,desig);

    // Format
    format_tle(orb,line1,line2);
    printf("%s\n%s\n%s\n",line0,line1,line2);
  }
  fclose(file);
  return 0;
}
