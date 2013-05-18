#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIM 81

int fgetline(FILE *,char *,int);
void rtrim(char *);

int main(int argc,char *argv[])
{
  char line[LIM];
  FILE *fitsfile;

  // Usage
  if (argc<2) {
    printf("Usage: %s <fitsfile>\n",argv[0]);
    printf("\n\nOutputs the header of <fitsfile>\n");
    return 1;
  }

  // Open file
  fitsfile=fopen(argv[1],"r");

  // Loop over file and output
  while (fgetline(fitsfile,line,LIM)>0) {
    rtrim(line);
    printf("%s\n",line);
    if (strcmp(line,"END")==0) break;
  }

  // Close file
  fclose(fitsfile);

  return 0;
}

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
}

// Removes trailing blanks from string s
void rtrim(char *s)
{
  int i,j=0,n;

  n=strlen(s);
  for (i=n;i>=0;i--)
    if (s[i]!='\0' && s[i]!='\n')
      if (!isspace(s[i]) && j==0) j=i;
  s[++j]='\0';

  return;
}
