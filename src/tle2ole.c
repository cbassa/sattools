#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#define LIM 256

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

int main(int argc,char *argv[])
{
  int reverse=0;
  char line[LIM],pline[LIM],tlefile[LIM];
  char line0[70],line1[70],line2[70];
  FILE *file;
  int arg=0;

  // Decode options
  while ((arg=getopt(argc,argv,"c:rh"))!=-1) {
    switch(arg) {

    case 'c':
      strcpy(tlefile,optarg);
      break;

    case 'r':
      reverse=1;
      break;

    case 'h':
    default:
      printf("-c <catalog>\n-r reverse\n");
      return 0;
    }
  }

  if (reverse==0) {
    file=fopen(tlefile,"r");
    while (fgetline(file,line,LIM)>0) {
      if (line[0]=='1') {
	strcpy(line0,pline);
	strcpy(line1,line);
	fgetline(file,line,LIM);
	strcpy(line2,line);
	printf("%s | %s | %s\n",line1,line2,line0);
      } else {
	strcpy(pline,line);
      }
    }
    fclose(file);
  } else {
    file=fopen(tlefile,"r");
    while (fgetline(file,line,LIM)>0) {
      printf("%.70s\n%.70s\n%.70s\n",line+144,line,line+72);
    }
    fclose(file);
  }

  return 0;
}
