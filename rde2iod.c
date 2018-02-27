#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define LIM 128

int fgetline(FILE *file,char *s,int lim);
int find_satno(char *desig0)
{
  FILE *file;
  int satno=99999,status;
  char desig[16];
  char *env,filename[LIM];

  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/desig.txt",env);
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"Designation file not found!\n");
    exit(0);
  }
  while (!feof(file)) {
    status=fscanf(file,"%d %s",&satno,desig);
    if (strcmp(desig,desig0)==0) 
      break;
  }
  fclose(file);

  return satno;
}


int main(int argc,char *argv[])
{
  FILE *file;
  char line[LIM];
  int intidy,intido,piece,site,year,month,day,hour,min,sec,fsec,satno;
  char desig[16],pdesig[16];
  int format,epoch,dummy,icsec;
  float csec,cang;
  int rah,ram,rafm,ded,dem,defm;
  float tm,tx,am,ax;
  char sign;
  int lineno=0;

  file=fopen(argv[1],"r");
  while (fgetline(file,line,LIM)>0) {
    if (strncmp(line,"2420",4)==0 && lineno==0) {
      lineno++;
      sscanf(line,"%04d %02d%02d",&site,&year,&month);
      sscanf(line+12,"%1d%1d%1d %3d%1d",&icsec,&dummy,&format,&dummy,&epoch);

      csec=0.1*icsec;
      cang=dummy;

      // Year switch
      if (year>50)
	year+=1900;
      else
	year+=2000;

      // Time accuracy
      tx=floor(log10(csec))+8;
      tm=floor(csec/pow(10.0,tx-8));

      // angle accuracy
      ax=floor(log10(cang))+8;
      am=floor(cang/pow(10.0,ax-8));

      if (ax>9.0) {
	ax=9.0;
	am=9.0;
      }

      continue;
    }
    if (strlen(line)<5 && lineno==1) {
      sscanf(line,"%d",&day);
      if (day==999)
	break;
      continue;
    }
    
    // Skip wrong lines
    if (!isdigit(line[0])) 
      continue;
    // Skip short lines
    if (strlen(line)<31)
      continue;

    // Scan line
    sscanf(line,"%02d%03d%02d",&intidy,&intido,&piece);
    sscanf(line+8,"%02d%02d%02d.%02d",&hour,&min,&sec,&fsec);
    sscanf(line+18,"%02d%02d%d",&rah,&ram,&rafm);
    sscanf(line+24,"%c%02d%02d%d",&sign,&ded,&dem,&defm);
    fsec*=10.0;
    // Format designation
    if (piece<26) {
      sprintf(desig,"%02d %03d%c",intidy,intido,piece+'A'-1);
      sprintf(pdesig,"%02d%03d%c",intidy,intido,piece+'A'-1);
    } else {
      fprintf(stderr,"Failed to understand designation!\n");
      fprintf(stderr,"%s\n",line);
      continue;
    }

    // Test data format
    if (format!=1) {
      fprintf(stderr,"Angle format %d not implemented!\n",format);
      fprintf(stderr,"%s\n",line);
      continue;
    }

    // Fractional RA
    if (rafm<10)
      rafm*=100;
    else if (rafm<100)
      rafm*=10;

    // Fractional DE
    if (defm<10)
      defm*=10;
    else if (defm<100)
      defm*=1;

    // Get satellite number
    satno=find_satno(pdesig);

    // Format IOD line
    printf("%05d %s   %04d G %04d%02d%02d%02d%02d%02d%03d %1.0f%1.0f %d%d ",satno,desig,site,year,month,day,hour,min,sec,fsec,tm,tx,format,epoch);
    printf("%02d%02d%03d%c%02d%02d%02d %1.0f%1.0f\n",rah,ram,rafm,sign,ded,dem,defm,am,ax);
  }
  fclose(file);

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

