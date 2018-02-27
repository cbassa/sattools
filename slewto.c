#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <netdb.h>
#include <sys/socket.h>

#include <netinet/in.h>
#include <netinet/in_systm.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>

#define LIM 128
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define PORT 7264
#define IP "127.0.0.1"

struct map {
  double alpha0,delta0,ra0,de0,azi0,alt0,q;
  char orientation[LIM];
  char nfd[LIM];
  char datadir[LIM],observer[32];
  double lat,lng;
  double mjd;
  float alt,timezone;
  int site_id;
} m;

void usage(void)
{

  return;
}

// Compute Date from Julian Day
void mjd2nfd(double mjd,char *nfd)
{
  double f,jd,dday;
  int z,alpha,a,b,c,d,e;
  int year,month,day,hour,min;
  float sec,x;

  jd=mjd+2400000.5;
  jd+=0.5;

  z=floor(jd);
  f=fmod(jd,1.);

  if (z<2299161)
    a=z;
  else {
    alpha=floor((z-1867216.25)/36524.25);
    a=z+1+alpha-floor(alpha/4.);
  }
  b=a+1524;
  c=floor((b-122.1)/365.25);
  d=floor(365.25*c);
  e=floor((b-d)/30.6001);

  dday=b-d-floor(30.6001*e)+f;
  if (e<14)
    month=e-1;
  else
    month=e-13;

  if (month>2)
    year=c-4716;
  else
    year=c-4715;

  day=(int) floor(dday);
  x=24.0*(dday-day);
  x=3600.*fabs(x);
  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  hour=x;
  sec=floor(1000.0*sec)/1000.0;

  sprintf(nfd,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);

  return;
}

// Send new position to telescope
void send_position(char *sra,char *sde)
{
  int skt;
  struct hostent *he;
  struct sockaddr_in addr;
  char packet[2048];
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
  //  file=fopen("/media/video/satobs/control/state.txt","w");
  //  if (file!=NULL) {
  //    fprintf(file,"restart");
  //    fclose(file);
  //  }

  // Set position
  //  file=fopen("/media/video/satobs/control/position.txt","w");
  //  if (file!=NULL) {
  //    fprintf(file,"%s %s\n",sra,sde);
  //    fclose(file);
  //  }

  return;
}

// Get observing site
void get_site(int site_id)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64],filename[LIM];

  sprintf(filename,"%s/data/sites.txt",m.datadir);
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with site information not found!\n");
    return;
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip
    if (strstr(line,"#")!=NULL)
      continue;

    // Strip newline
    line[strlen(line)-1]='\0';

    // Read data
    sscanf(line,"%4d %2s %lf %lf %f",
	   &id,abbrev,&lat,&lng,&alt);
    strcpy(observer,line+38);

    // Change to km
    alt/=1000.0;
    
    if (id==site_id) {
      m.lat=lat;
      m.lng=lng;
      m.alt=alt;
      m.site_id=id;
      strcpy(m.observer,observer);
    }

  }
  fclose(file);

  return;
}

void initialize(void)
{
  int i;
  char *env;

  // Default Map parameters
  m.azi0=0;
  m.alt0=90.0;

  m.lat=0.0;
  m.lng=0.0;
  m.alt=0.0;

  // Default settings
  m.site_id=0;

  // Get environment variables
  env=getenv("ST_DATADIR");
  if (env!=NULL) {
    strcpy(m.datadir,env);
  } else {
    printf("ST_DATADIR environment variable not found.\n");
  }
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    get_site(atoi(env));
  } else {
    printf("ST_COSPAR environment variable not found.\n");
  }

  return;
}


// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}


// Convert Sexagesimal into Decimal
double sex2dec(char *s)
{
  double x;
  float deg,min,sec;
  char t[LIM];

  strcpy(t,s);

  deg=fabs(atof(strtok(t," :")));
  min=fabs(atof(strtok(NULL," :")));
  sec=fabs(atof(strtok(NULL," :")));

  x=(double) deg+(double) min/60.+(double) sec/3600.;
  if (s[0]=='-') x= -x;

  return x;
}

// Compute Julian Day from Date
double date2mjd(int year,int month,double day)
{
  int a,b;
  double jd;

  if (month<3) {
    year--;
    month+=12;
  }

  a=floor(year/100.);
  b=2.-a+floor(a/4.);

  if (year<1582) b=0;
  if (year==1582 && month<10) b=0;
  if (year==1582 && month==10 && day<=4) b=0;

  jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

  return jd-2400000.5;
}

// Present nfd
void nfd_now(char *s)
{
  time_t rawtime;
  struct tm *ptm;

  // Get UTC time
  time(&rawtime);
  ptm=gmtime(&rawtime);
    
  sprintf(s,"%04d-%02d-%02dT%02d:%02d:%02d",ptm->tm_year+1900,ptm->tm_mon+1,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
  
  return;
}

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min,sec;
  double mjd,dday;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%02d",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;

  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Convert horizontal into equatorial coordinates
void horizontal2equatorial(double mjd,double azi,double alt,double *ra,double *de)
{
  double h;

  h=atan2(sin(azi*D2R),cos(azi*D2R)*sin(m.lat*D2R)+tan(alt*D2R)*cos(m.lat*D2R))*R2D;
  *ra=modulo(gmst(mjd)+m.lng-h,360.0);
  *de=asin(sin(m.lat*D2R)*sin(alt*D2R)-cos(m.lat*D2R)*cos(alt*D2R)*cos(azi*D2R))*R2D;
  if (*ra<0.0)
    *ra+=360.0;

  return;
}

// Convert equatorial into horizontal coordinates
void equatorial2horizontal(double mjd,double ra,double de,double *azi,double *alt)
{
  double h;

  h=gmst(mjd)+m.lng-ra;
  
  *azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(m.lat*D2R)-tan(de*D2R)*cos(m.lat*D2R))*R2D,360.0);
  *alt=asin(sin(m.lat*D2R)*sin(de*D2R)+cos(m.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

  return;
}

double parallactic_angle(double mjd,double ra,double de)
{
  double h,q;

  h=gmst(mjd)+m.lng-ra;

  q=atan2(sin(h*D2R),(tan(m.lat*D2R)*cos(de*D2R)-sin(de*D2R)*cos(h*D2R)))*R2D;

  return q;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int f,int len)
{
  int i;
  double sec,deg,min;
  char sign;
  char *form[]={"::",",,","hms","  "};

  sign=(x<0 ? '-' : ' ');
  x=3600.*fabs(x);

  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;

  if (len==7) sprintf(s,"%c%02i%c%02i%c%07.4f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==6) sprintf(s,"%c%02i%c%02i%c%06.3f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==5) sprintf(s,"%c%02i%c%02i%c%05.2f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==4) sprintf(s,"%c%02i%c%02i%c%04.1f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==2) sprintf(s,"%c%02i%c%02i%c%02i%c",sign,(int) deg,form[f][0],(int) min,form[f][1],(int) floor(sec),form[f][2]);

  return;
}


int main(int argc,char *argv[])
{
  int arg=0,haflag=0,dry=0;
  char sra[16],sde[16];
  double ha;
  FILE *file;

  // Initialize
  initialize();

  // Get time
  nfd_now(m.nfd);
  m.mjd=nfd2mjd(m.nfd);

    // Decode options
  while ((arg=getopt(argc,argv,"m:t:H:R:D:A:E:hn"))!=-1) {
    switch(arg) {

    case 'n': 
      dry=1;
      break;

    case 't':
      strcpy(m.nfd,optarg);
      m.mjd=nfd2mjd(m.nfd);
      break;

    case 'm':
      m.mjd=atof(optarg);
      mjd2nfd(m.mjd,m.nfd);
      break;

    case 'H':
      ha=atof(optarg);
      haflag=1;
      strcpy(m.orientation,"equatorial");
      break;

    case 'R':
      m.ra0=15.0*sex2dec(optarg);
      strcpy(m.orientation,"equatorial");
      break;

    case 'D':
      m.de0=sex2dec(optarg);
      strcpy(m.orientation,"equatorial");
      break;

    case 'A':
      m.azi0=modulo(atof(optarg)+180.0,360.0);
      strcpy(m.orientation,"horizontal");
      break;

    case 'E':
      m.alt0=atof(optarg);
      strcpy(m.orientation,"horizontal");
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

  // Compute RA from HA
  if (haflag==1) {
    m.ra0=modulo(gmst(m.mjd)+m.lng-ha,360.0);
  } else {
    ha=modulo(gmst(m.mjd)+m.lng-m.ra0,360.0);
    if (ha>180.0)
      ha-=360.0;
  }
  // Compute RA and Dec
  if (strcmp(m.orientation,"horizontal")==0) 
    horizontal2equatorial(m.mjd,m.azi0,m.alt0,&m.ra0,&m.de0);
  else if (strcmp(m.orientation,"equatorial")==0) 
    equatorial2horizontal(m.mjd,m.ra0,m.de0,&m.azi0,&m.alt0);

  // Parallactic angle
  m.q=parallactic_angle(m.mjd,m.ra0,m.de0);

  // Get sexagesimal
  dec2sex(m.ra0/15.0,sra,0,5);
  dec2sex(m.de0,sde,0,4);

  // Print to screen
  printf("%s R:%s D: %s H: %7.3f A: %6.3f E: %6.3f q: %5.2f\n",m.nfd,sra,sde,ha,modulo(m.azi0-180.0,360.0),m.alt0,m.q);

  if (dry==0) {
    // Send position
    send_position(sra,sde);
    
    // Log
    file=fopen("position.txt","a");
    fprintf(file,"%s %lf %lf %f\n",m.nfd,m.ra0,m.de0,m.q);
    fclose(file);
  }

  return 0;
}
