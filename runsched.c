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

#define PORT 7624
#define MAXPORT 7674
#define IP "127.0.0.1"
#define LIM 2048
#define NMAX 256

#define SCHEDULED 0
#define STARTED 1
#define FINISHED 2

// Observation struct contains observation time, celestial coords and camera name
struct observation {
  char stime[32],sra[32],sde[32],camname[15],startstop[10];
  time_t ptime;
  float dt;
};

int fgetline(FILE *file,char *s,int lim);
void send_position(char *sra,char *sde,char *datadir,char *obsdir,char *camname);
void stop_obs(char *datadir,char *obsdir,char *camname);
time_t decode_time(char *stm);

int main(int argc, char *argv[]) 
{
  int i=0,nobs,flag=0;
  time_t rawtime,aimtime;
  struct tm *ptm,*rtm;
  char buf[20],line[LIM],stm[20],sra[32],sde[32],pra[32],pde[32];
  FILE *file;
  struct observation obs[NMAX];
  char *env;
  char datadir[128],obsdir[128];
	int nextobs, dtnext;

  // Get environment variables
  env=getenv("ST_DATADIR");
  if (env!=NULL) {
    strcpy(datadir,env);
  } else {
    printf("ST_DATADIR environment variable not found.\n");
  }

  // Get environment variables
  env=getenv("ST_OBSDIR");
  if (env!=NULL) {
    strcpy(obsdir,env);
  } else {
    printf("ST_OBSDIR environment variable not found.\n");
  }

  // For ever loop
  for (;;) {
    // Read file
    i=0;
    file=fopen("schedule.txt","r");
    while (fgetline(file,line,LIM)>0) {
      sscanf(line,"%s %s %s %s %s",obs[i].stime,obs[i].sra,obs[i].sde,obs[i].camname,obs[i].startstop);
      obs[i].ptime=decode_time(obs[i].stime);
      
      i++;
    }
    fclose(file);
    nobs=i;

    // Get local time
    time(&rawtime);
    //rawtime-=3600;
    

    // Print UTC time
    ptm=gmtime(&rawtime);
    strftime(buf,32,"%Y-%m-%dT%H:%M:%S",ptm);

    // Make raw time UTC to compare with scheduled time
    rawtime=mktime(ptm);

    // Show current raw time, just to check
    // printf("%s\n",ctime(&rawtime));

    // Compute time differences
    for (i=0;i<nobs;i++) {
      obs[i].dt=difftime(obs[i].ptime,rawtime);
    }

    nextobs=-1;
    dtnext=9999999;
    // Loop over observations
    for (i=0;i<nobs;i++) {
      if (obs[i].dt>0.0) {
	if(obs[i].dt < dtnext){
	  nextobs=i;
	  dtnext=obs[i].dt;
	}
	//				printf("%4.0f %s %s %s %s\n",obs[i].dt,obs[i].stime,obs[i].sra,obs[i].sde,obs[i].startstop);
	//				break;
      } else if (obs[i].dt==0) {
	if(strstr(obs[i].startstop,"tart")!=NULL){
	  //printf("Slewing to %s %s\n",obs[i].sra,obs[i].sde);
	  send_position(obs[i].sra,obs[i].sde,datadir,obsdir,obs[i].camname);
	} else if(strstr(obs[i].startstop,"top")!=NULL){
	  stop_obs(datadir,obsdir,obs[i].camname);
	}
      }
    }
    
    if(nextobs>=0){
      // print next observation data if any found
      printf("%4.0f %s %s %s %s\n",obs[nextobs].dt,obs[nextobs].stime,obs[nextobs].sra,obs[nextobs].sde,obs[nextobs].startstop);
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

// Read data/cameras.txt in search of specified camera name and return complete camera details line
int read_cameras(char *camname,char *datadir,char *camera)
{
  FILE *file;
  char line[127],filename[127];

  sprintf(filename,"%s/data/cameras.txt",datadir);
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with camera information not found!\n");
    return -1;
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip commented lines
    if (strstr(line,"#")!=NULL)
      continue;
    if(strstr(line, camname)!=NULL){
      strcpy(camera, line);
      return 0;
    }
  }
  fclose(file);
  return -1;
}


// Send new position to telescope
void send_position(char *sra,char *sde,char *datadir,char *obsdir,char *camname)
{
  int skt, port;
  struct hostent *he;
  struct sockaddr_in addr;
  char packet[LIM];
  FILE *file;
  float ra,de;
  char camera[128],fname[128];
  float f;
  char s[31];

  // Check if camera is fixed
  // read complete line from data/cameras.txt describing the scheduled camera
  read_cameras(camname,datadir,camera);  // search for camera name
  sscanf(camera,"%s %f %f %f %f %s", s, &f, &f, &f, &f, s);
  // Look for "fix" string to jump over slewing routines.
  if(strstr(s,"ix")==NULL){
  
    // Old packet style
    //  sprintf(packet,"<newNumberVector device='Celestron GPS' name='EQUATORIAL_EOD_COORD_REQUEST'><oneNumber name='RA'>%s</oneNumber><oneNumber name='DEC'>%s</oneNumber></newNumberVector>",sra,sde);

    // New packet style (as of 2013-08-20)
    sprintf(packet,"<newNumberVector device='Celestron GPS' name='EQUATORIAL_EOD_COORD'><oneNumber name='RA'>%s</oneNumber><oneNumber name='DEC'>%s</oneNumber></newNumberVector>",sra,sde);

		printf("Slewing to %s %s\n",sra,sde);

    // Send TCP packet
    skt=socket(AF_INET,SOCK_STREAM,0);
    addr.sin_family=AF_INET;
    port=PORT;
    addr.sin_port=htons(port);
    he=gethostbyname(IP);
    bcopy(he->h_addr,(struct in_addr *) &addr.sin_addr,he->h_length);
    while((connect(skt,(struct sockaddr *) &addr,sizeof(addr))<0) && (port < MAXPORT)) {
      fprintf(stderr,"Connection refused by remote host on port %04d.\n",port);
      port++;
      // Skip port 7265... used by some other service?
      if(port==7265) port++;
      fprintf(stderr,"Trying port %04d.\n",port);

      addr.sin_port=htons(port);
      he=gethostbyname(IP);
      bcopy(he->h_addr,(struct in_addr *) &addr.sin_addr,he->h_length);

    }
    if(port>=MAXPORT) return;

    printf("Connected to Indi server on port %04d.\n",port);

    write(skt,packet,strlen(packet));
    close(skt); 
  }

	printf("Starting new observation\n");
   
  // Set restart
  sprintf(fname,"%s/control/state.txt",obsdir);
  file=fopen(fname,"w");
  if (file!=NULL) {
    fprintf(file,"restart");
    fclose(file);
  }

  // Set position
  sprintf(fname,"%s/control/position.txt",obsdir);
  file=fopen(fname,"w");
  if (file!=NULL) {
    fprintf(file,"%s %s\n",sra,sde);
    fclose(file);
  }

  // Set camera
  sprintf(fname,"%s/control/camera.txt",obsdir);
  file=fopen(fname,"w");
  if (file!=NULL) {
    fprintf(file,"%s",camera);
    fclose(file);
  }

  return;
}


// Send stop observation signal
void stop_obs(char *datadir,char *obsdir,char *camname)
{
  FILE *file;
  char camera[128],fname[128];
  float f;
  char s[31];

  // Retrieve Camera data
  // read complete line from data/cameras.txt describing the scheduled camera
  read_cameras(camname,datadir,camera);  // search for camera name
  sscanf(camera,"%s %f %f %f %f %s", s, &f, &f, &f, &f, s);

	printf("Stop observation\n");
   
  // Set stop
  sprintf(fname,"%s/control/state.txt",obsdir);
  file=fopen(fname,"w");
  if (file!=NULL) {
    fprintf(file,"stop");
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
