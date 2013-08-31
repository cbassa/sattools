#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <netdb.h>
#include <sys/socket.h>

#include <netinet/in.h>
#include <netinet/in_systm.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>

#define PORT 7264
#define IP "127.0.0.1"
#define LIM 2048
#define NMAX 128

#define SCHEDULED 0
#define STARTED 1
#define FINISHED 2

struct observation {
  char stime[20],sra[15],sde[15];
  time_t ptime;
  float dt;
};

int fgetline(FILE *file,char *s,int lim);
void send_position(char *sra,char *sde);
time_t decode_time(char *stm);


int main(int argc, char *argv[]) 
{
  int i=0,nobs,flag=0;
  time_t rawtime,aimtime;
  struct tm *ptm,*rtm;
  char buf[20],line[LIM],stm[20],sra[15],sde[15],pra[15],pde[15];
  FILE *file;
  struct observation obs[NMAX];

  // For ever loop
  for (;;) {
    // Read file
    i=0;
    file=fopen("schedule.txt","r");
    while (fgetline(file,line,LIM)>0) {
      sscanf(line,"%s %s %s",obs[i].stime,obs[i].sra,obs[i].sde);
      obs[i].ptime=decode_time(obs[i].stime);
      
      i++;
    }
    fclose(file);
    nobs=i;

    // Get local time
    time(&rawtime);
    rawtime-=3600;

    // Print UTC time
    ptm=gmtime(&rawtime);
    strftime(buf,20,"%Y-%m-%dT%H:%M:%S",ptm);

    // Compute time differences
    for (i=0;i<nobs;i++) 
      obs[i].dt=difftime(obs[i].ptime,rawtime);

    // Loop over observations
    for (i=0;i<nobs;i++) {
      if (obs[i].dt>0.0) {
	printf("%4.0f %s %s %s\n",obs[i].dt,obs[i].stime,obs[i].sra,obs[i].sde);
	break;
      } else if (obs[i].dt==0) {
	printf("Slewing to %s %s\n",obs[i].sra,obs[i].sde);
	send_position(obs[i].sra,obs[i].sde);
      }
    }

    // Sleep
    sleep(1);
  }

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

// Send new position to telescope
void send_position(char *sra,char *sde)
{
  int skt;
  struct hostent *he;
  struct sockaddr_in addr;
  char packet[LIM];
  FILE *file;
  float ra,de;


  // Old packet style
  //  sprintf(packet,"<newNumberVector device='Celestron GPS' name='EQUATORIAL_EOD_COORD_REQUEST'><oneNumber name='RA'>%s</oneNumber><oneNumber name='DEC'>%s</oneNumber></newNumberVector>",sra,sde);

  // New packet style (as of 2013-08-20)
  sprintf(packet,"<newNumberVector device='Celestron GPS' name='EQUATORIAL_EOD_COORD'><oneNumber name='RA'>%s</oneNumber><oneNumber name='DEC'>%s</oneNumber></newNumberVector>",sra,sde);

  // Send TCP packet
  skt=socket(AF_INET,SOCK_STREAM,0);
  addr.sin_family=AF_INET;
  addr.sin_port=htons(PORT);
  he=gethostbyname(IP);
  bcopy(he->h_addr,(struct in_addr *) &addr.sin_addr,he->h_length);
  if(connect(skt,(struct sockaddr *) &addr,sizeof(addr))<0) {
    fprintf(stderr,"Connection refused by remote host.\n");
    return;
  }
  write(skt,packet,strlen(packet));
  close(skt); 
 
  // Set restart
  file=fopen("/media/video/satobs/control/state.txt","w");
  if (file!=NULL) {
    fprintf(file,"restart");
    fclose(file);
  }

  // Set position
  file=fopen("/media/video/satobs/control/position.txt","w");
  if (file!=NULL) {
    fprintf(file,"%s %s\n",sra,sde);
    fclose(file);
  }

  return;
}

// Decode time
time_t decode_time(char *stm)
{
  time_t aimtime;
  struct tm *rtm;
  int d;

  rtm=gmtime(&aimtime);
  sscanf(stm,"%04d-%02d-%02dT%02d:%02d:%02d",&rtm->tm_year,&rtm->tm_mon,&rtm->tm_mday,&rtm->tm_hour,&rtm->tm_min,&rtm->tm_sec);
  rtm->tm_year-=1900;
  rtm->tm_mon--;
  aimtime=mktime(rtm);

  return aimtime;
}
