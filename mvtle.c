#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
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

// Format TLE
void format_tle(orbit_t orb,char *line1,char *line2)
{
  int i,csum;
  char sbstar[]=" 00000-0",bstar[13];

  // Format Bstar term
  if (fabs(orb.bstar)>1e-9) {
    sprintf(bstar,"%11.4e",10*orb.bstar);
    sbstar[0] = bstar[0];  sbstar[1] = bstar[1];  sbstar[2] = bstar[3];  sbstar[3] = bstar[4];
    sbstar[4] = bstar[5];  sbstar[5] = bstar[6];  sbstar[6] = bstar[8];  sbstar[7] = bstar[10];  sbstar[8] = '\0';
  }
  // Print lines
  sprintf(line1,"1 %05dU %-8s %2d%012.8f  .00000000  00000-0 %8s 0    0",orb.satno,orb.desig,orb.ep_year-2000,orb.ep_day,sbstar);
  sprintf(line2,"2 %05d %8.4f %8.4f %07.0f %8.4f %8.4f %11.8f%5ld",orb.satno,DEG(orb.eqinc),DEG(orb.ascn),1E7*orb.ecc,DEG(orb.argp),DEG(orb.mnan),orb.rev,orb.norb);

  // Compute checksums
  for (i=0,csum=0;i<strlen(line1);i++) {
    if (isdigit(line1[i]))
      csum+=line1[i]-'0';
    else if (line1[i]=='-')
      csum++;
  }
  sprintf(line1,"%s%d",line1,csum%10);
  for (i=0,csum=0;i<strlen(line2);i++) {
    if (isdigit(line2[i]))
      csum+=line2[i]-'0';
    else if (line2[i]=='-')
      csum++;
  }
  sprintf(line2,"%s%d",line2,csum%10);

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
  while ((arg=getopt(argc,argv,"c:i:I:D:N:"))!=-1) {
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
      return 0;
      break;

    default:
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
