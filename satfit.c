#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "cel.h"
#include "sgdp4h.h"
#include <getopt.h>

#define LIM 80
#define NMAX 256
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

struct point {
  int flag,satno;
  char desig[10];
  double mjd,ra,de,rac,dec;
  float st,sr;
  char iod_line[LIM];
  double dx,dy,dr,dt;
  xyz_t obspos;
};
struct data {
  int n,nsel;
  struct point *p;
  double chisq,rms;
} d;
struct site {
  int id;
  double lng,lat;
  float alt;
  char observer[64];
};
orbit_t orb;
struct site get_site(int site_id);
struct point decode_iod_observation(char *iod_line);
int fgetline(FILE *file,char *s,int lim);
double modulo(double x,double y);
double gmst(double mjd);
double dgmst(double mjd);
double date2mjd(int year,int month,double day);
double mjd2doy(double mjd,int *yr);
void mjd2date(double mjd,int *year,int *month,double *day);
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel);
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de);
void forward(double ra0,double de0,double ra,double de,double *x,double *y);
struct data read_data(char *filename);
void versafit(int m,int n,double *a,double *da,double (*func)(double *),double dchisq,double tol,char *opt);
double chisq(double a[]);
orbit_t read_tle(char *filename,int satno);
void format_tle(orbit_t orb,char *line1,char *line2);
void highlight(float x0,float y0,float x,float y,int flag);
void time_range(double *mjdmin,double *mjdmax,int flag);
void print_tle(orbit_t orb,char *filename);
void fit(orbit_t orb,int *ia);
void usage();

xyz_t get_position(double r0,int i0)
{
  int i;
  double rr,drr,r;
  xyz_t pos;
  double x,y,z;

  // Initial range
  rr=100.0;
  do {
    x=d.p[i0].obspos.x+rr*cos(d.p[i0].ra*D2R)*cos(d.p[i0].de*D2R);
    y=d.p[i0].obspos.y+rr*sin(d.p[i0].ra*D2R)*cos(d.p[i0].de*D2R);
    z=d.p[i0].obspos.z+rr*sin(d.p[i0].de*D2R);
    r=sqrt(x*x+y*y+z*z);
    drr=r-r0;
    rr-=drr;
  } while (fabs(drr)>0.01);
  
  pos.x=x;
  pos.y=y;
  pos.z=z;

  return pos;
}

int period_search(void)
{
  int i,j,i1,i2;
  float dt;
  int nrev,nrevmin,nrevmax;
  char line1[70],line2[70];
  int ia[7];

  // Set fitting parameters
  for (i=0;i<6;i++)
    ia[i]=1;
  ia[6]=0;

  // Select all points
  for (i=0;i<d.n;i++)
    d.p[i].flag=2;
  

  // Print observations
  printf("Present observations:\n");
  for (i=0;i<d.n;i++) 
    printf("%3d: %s\n",i+1,d.p[i].iod_line);
  printf("\nSelect center observations of both arcs: ");
  scanf("%d %d",&i1,&i2);
  dt=d.p[i2].mjd-d.p[i1].mjd;
  printf("\nTime passed: %f days\n",dt);
  printf("Provide revolution range to search over [min,max]: ");
  scanf("%d %d",&nrevmin,&nrevmax);

  for (nrev=nrevmin;nrev<nrevmax+1;nrev++) {
    orb.satno=79000+nrev;
    orb.rev=nrev/dt;
    
    // Set parameters
    for (i=0;i<7;i++) 
      ia[i]=0;
	
    // Loop over parameters
    for (i=0;i<6;i++) {
      if (i==0) ia[4]=1;
      if (i==1) ia[1]=1;
      if (i==2) ia[0]=1;
      if (i==3) ia[3]=1;
      if (i==4) ia[2]=1;
      if (i==5) ia[5]=1;

      for (j=0;j<5;j++)
	fit(orb,ia);
    }

    format_tle(orb,line1,line2);
    printf("%s\n%s\n# %d revs, %f revs/day, %f\n",line1,line2,nrev,nrev/dt,d.rms);
  }


  return;
}

int psearch(void)
{
  int i,satno=99300;
  double mjdmin,mjdmax;
  int ia[7]={0,0,0,0,0,0,0};
  double ecc,eccmin,eccmax,decc;
  double rev,revmin,revmax,drev;
  double argp,argpmin,argpmax,dargp;
  char line1[70],line2[70];
  FILE *file;

  // Provide 
  printf("Mean motion [min, max, stepsize]: \n");
  scanf("%lf %lf %lf",&revmin,&revmax,&drev);
  printf("Eccentricity [min, max, stepsize]: \n");
  scanf("%lf %lf %lf",&eccmin,&eccmax,&decc);
  //  printf("Argument of perigee [min, max, stepsize]: \n");
  //  scanf("%lf %lf %lf",&argpmin,&argpmax,&dargp);


  // Step 1: select all points
  //  for (i=0;i<d.n;i++)
  //    d.p[i].flag=2;

  // Step 2: get time range
  time_range(&mjdmin,&mjdmax,2);

  file=fopen("search.dat","w");

  // Step 4: Loop over eccentricity
  for (rev=revmin;rev<revmax;rev+=drev) {
    for (ecc=eccmin;ecc<eccmax;ecc+=decc) {
      //      for (argp=argpmin;argp<argpmax;argp+=dargp) {
	orb.satno=satno;
	orb.ecc=ecc;
	orb.rev=rev;
	//orb.argp=argp*D2R;
    
	// Set parameters
	for (i=0;i<7;i++) 
	  ia[i]=0;
	
	// Step 4: loop over parameters
	for (i=0;i<5;i++) {
	  if (i==0) ia[4]=1;
	  if (i==1) ia[1]=1;
	  if (i==2) ia[0]=1;
	  //	if (i==3) ia[5]=1;
	  if (i==4) ia[3]=1;
	  
	  // Do fit
	  fit(orb,ia);
	}
	fit(orb,ia);
	fit(orb,ia);
	fit(orb,ia);
	printf("%8.5lf %8.6lf %8.3lf %8.3lf %8.3lf %8.3lf %8.5lf\n",orb.rev,orb.ecc,orb.argp*R2D,orb.ascn*R2D,orb.mnan*R2D,orb.eqinc*R2D,d.rms);
	fprintf(file,"%8.5lf %8.6lf %8.3lf %8.3lf %8.3lf %8.3lf %8.5lf\n",orb.rev,orb.ecc,orb.argp*R2D,orb.ascn*R2D,orb.mnan*R2D,orb.eqinc*R2D,d.rms);
      }
      fprintf(file,"\n");
      //    }
  }
  fclose(file);

  return orb.satno;
}


int circular_fit(void)
{
  int i;
  double mjdmin,mjdmax;
  int ia[7]={0,0,0,0,0,0,0};

  // Step 1: select all points
  //  for (i=0;i<d.n;i++)
  //    d.p[i].flag=2;

  // Step 2: get time range
  time_range(&mjdmin,&mjdmax,2);

  // Step 3: set initial orbit
  orb.satno=d.p[0].satno;
  orb.eqinc=0.5*M_PI;
  orb.ascn=0.0;
  orb.ecc=0.0;
  orb.argp=0.0;
  orb.mnan=0.0;
  orb.rev=14.0;
  orb.bstar=0.5e-4;
  orb.ep_day=mjd2doy(0.5*(mjdmin+mjdmax),&orb.ep_year);

  // Step 4: loop over parameters
  for (i=0;i<4;i++) {
    if (i==0) ia[4]=1;
    if (i==1) ia[1]=1;
    if (i==2) ia[0]=1;
    if (i==3) ia[5]=1;

    // Do fit
    fit(orb,ia);
  }
  fit(orb,ia);

  return orb.satno;
}

int adjust_fit(void)
{
  int i;
  double mjdmin,mjdmax;
  int ia[6]={0,0,0,0,0,0};

  // Step 1: select all points
  for (i=0;i<d.n;i++)
    d.p[i].flag=2;

  // Step 2: loop over parameters
  for (i=0;i<2;i++) {
    if (i==0) ia[4]=1;
    if (i==1) ia[1]=1;

    // Do fit
    fit(orb,ia);
  }
  fit(orb,ia);

  return orb.satno;
}

void old_circular_fit(void)
{
  int i,j,i0,i1;
  float r0=6500;
  xyz_t pos0,pos1;
  double ang,dt,w,w0;

  // Get end points
  for (i=0,j=0;i<d.n;i++) {
    if (d.p[i].flag==2) {
      if (j==0)
	i0=i;
      i1=i;
      j++;
    }
  }

  // Time difference
  dt=86400.0*(d.p[i1].mjd-d.p[i0].mjd);
  i=0;
  do {
    w0=360.0/(2.0*M_PI*sqrt(r0*r0*r0/398600));
    // Get positions
    pos0=get_position(r0,i0);
    pos1=get_position(r0,i1);
    
    // Compute angle
    ang=acos((pos0.x*pos1.x+pos0.y*pos1.y+pos0.z*pos1.z)/(r0*r0))*R2D;
    
    // Angular motion (deg/sec);
    w=ang/dt;

    r0+=1000.0*(w0-w);
    i++;
  } while (fabs(w0-w)>1e-5 && i<1000);
  printf("%f\n",r0);
  return;
}

int main(int argc,char *argv[])
{
  int i,j,nobs=0;
  int redraw=1,plot_residuals=0,adjust=0,quit=0;
  int ia[]={0,0,0,0,0,0,0};
  float dx[]={0.1,0.1,0.35,0.35,0.6,0.6,0.85},dy[]={0.0,-0.25,0.0,-0.25,0.0,-0.25,0.0};
  char c;
  int mode=0,posn=0,click=0;
  float x0,y0,x,y;
  float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
  char string[64],bstar[10]=" 50000-4",line0[72],line1[72],line2[72],text[10];
  char filename[64];
  int satno=-1;
  double mjdmin,mjdmax;
  int arg=0,elset=0,circular=0,tleout=0,noplot=0;
  char *datafile,*catalog,tlefile[LIM];
  orbit_t orb0;

  // Decode options
  while ((arg=getopt(argc,argv,"d:c:i:haCo:p"))!=-1) {
    switch(arg) {
    case 'd':
      datafile=optarg;
      break;

    case 'c':
      catalog=optarg;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'C':
      circular=1;
      break;

    case 'p':
      noplot=1;
      break;

    case 'o':
      tleout=1;
      strcpy(tlefile,optarg);
      break;

    case 'h':
      usage();
      return 0;
      break;

    case 'a':
      adjust=1;
      break;

    default:
      usage();
      return 0;
    }
  }

  // Read data
  d=read_data(datafile);
  time_range(&mjdmin,&mjdmax,1);

  // Read TLE
  if (satno>=0) {
    orb=read_tle(catalog,satno);
  }

  freopen("/tmp/stderr.txt","w",stderr);

  // Fit circular orbit
  if (circular==1) {
    for (i=0;i<d.n;i++)
      d.p[i].flag=2;
    satno=circular_fit();
    plot_residuals=1;
    quit=1;

    // Dump tle
    if (tleout==1) 
      print_tle(orb,tlefile);
  }

  // Adjust
  if (adjust==1) {
    orb0=orb;
    adjust_fit();
    fit(orb,ia);
    printf("%05d %8.3f %8.3f %8.3f %s %8.3f\n",satno,DEG(orb.mnan-orb0.mnan),DEG(orb.ascn-orb0.ascn),d.rms,datafile,mjdmin-(SGDP4_jd0-2400000.5));
    plot_residuals=1;
    redraw=1;
    quit=1;

    // Dump tle
    if (tleout==1) 
      print_tle(orb,tlefile);
  }

  // Exit before plotting
  if (quit==1 && noplot==1) {
    free(d.p);
    fclose(stderr);
    return 0;
  }

  cpgopen("/xs");
  cpgask(0);

  // For ever loop
  for (;;) {
    if (redraw==1) {
      cpgpage();
      cpgsvp(0.1,0.95,0.0,0.18);
      cpgswin(0.0,1.0,-0.5,0.5);

      // Buttons
      cpgtext(0.12,-0.05,"Inclination");
      cpgtext(0.372,-0.05,"Eccentricity");
      cpgtext(0.62,-0.05,"Mean Anomaly");
      cpgtext(0.87,-0.05,"B\\u*\\d");
      cpgtext(0.12,-0.3,"Ascending Node");
      cpgtext(0.37,-0.3,"Arg. of Perigee");
      cpgtext(0.62,-0.3,"Mean Motion");

      // Toggles
      for (i=0;i<7;i++) {
	cpgpt1(dx[i],dy[i],19);
	if (ia[i]==1) {
	  cpgsci(2);
	  cpgpt1(dx[i],dy[i],16);
	  cpgsci(1);
	}
      }

      // Plot map
      cpgsvp(0.1,0.9,0.2,0.9);
      cpgswin(xmax,xmin,ymin,ymax);
      cpgbox("BCTSN",0.,0,"BCTSN",0.,0);
      cpglab("Right Ascension","Declination"," ");

      if (satno>0) {
	// Plot tle
	format_tle(orb,line1,line2);
	cpgmtxt("T",2.0,0.0,0.0,line1);
	cpgmtxt("T",1.0,0.0,0.0,line2);
      }

      // Plot points
      for (i=0;i<d.n;i++) {
	if (d.p[i].flag>=1) {
	  cpgpt1(d.p[i].ra,d.p[i].de,17);
	  sprintf(text," %d",i+1);
	  cpgtext(d.p[i].ra,d.p[i].de,text);
	  if (plot_residuals==1) {
	    cpgmove(d.p[i].ra,d.p[i].de);
	    cpgdraw(d.p[i].rac,d.p[i].dec);
	  }
	  if (d.p[i].flag==2) {
	    cpgsci(2);
	    cpgpt1(d.p[i].ra,d.p[i].de,4);
	    cpgsci(1);
	  }
	}
      }
    }

    // Quit
    if (quit==1)
      break;

    // Get cursor
    cpgband(mode,posn,x0,y0,&x,&y,&c);

    // Help
    if (c=='h' || c=='?') {
      printf("q   Quit\nw   Write TLE file\nP   Search mean motion space\nf   Fit orbit\ns   Select observations in current window\nu   Unselect observations\n1-7 Toggle parameter\nC   Fit circular orbit\nS   Search mean motion/eccentricity space\nc   Change a parameter\nz   Zoom\nr   Unzoom\nt   Toggle starting orbits (LEO, GTO, GEO, HEO)\n");
    }

    // Quit
    if (c=='q' || c=='Q')
      break;
   
    // Period search
    if (c=='P') {
      period_search();
    }

    // Fit
    if (c=='f') {
      // Count points
      for (i=0,nobs=0;i<d.n;i++)
	if (d.p[i].flag==2)
	  nobs++;
      if (satno<0) {
	printf("No elements loaded!\n");
      } else if (nobs==0) {
	printf("No points selected!\n");
      } else {
	fit(orb,ia);
	printf("%d %.5f\n",nobs,d.rms);
	plot_residuals=1;
	redraw=1;
	continue;
      }
    }

    // Write TLE
    if (c=='w') {
      printf("TLE filename to write: ");
      scanf("%s",filename);
      print_tle(orb,filename);
      printf("\n================================================================================\n");
      continue;
    }

    // Highlight
    if (c=='s') {
      highlight(xmin,ymin,xmax,ymax,2);
      time_range(&mjdmin,&mjdmax,2);
      for (i=0,nobs=0;i<d.n;i++)
	if (d.p[i].flag==2)
	  nobs++;
      click=0;
      mode=0;
      redraw=1;
      continue;
    }

    // Unselect
    if (c=='U') {
      for (i=0;i<d.n;i++)
	d.p[i].flag=1;
      time_range(&mjdmin,&mjdmax,1);
      redraw=1;
      continue;
    }

    // Unselect
    if (c=='u') {
      for (i=0;i<d.n;i++) 
	if (d.p[i].flag==2)
	  d.p[i].flag=1;
      redraw=1;
      continue;
    }

    // Toggles
    if (isdigit(c) && c-'0'>=1 && c-'0'<8) {
      if (ia[c-49]==0) 
	ia[c-49]=1;
      else if (ia[c-49]==1) 
	ia[c-49]=0;
      redraw=1;
      continue;
    }

    // Circular fit
    if (c=='C') {
      satno=circular_fit();
      plot_residuals=1;
      printf("%.3f\n",d.rms);
      ia[0]=ia[1]=ia[4]=ia[5]=1;
      redraw=1;
    }

    // Search
    if (c=='S') {
      satno=psearch();
      plot_residuals=1;
      ia[0]=ia[1]=ia[4]=ia[5]=1;
      redraw=1;
    }

    // Change
    if (c=='c') {
      printf("(1) Inclination,     (2) Ascending Node,   (3) Eccentricity,\n(4) Arg. of Perigee, (5) Mean Anomaly,     (6) Mean Motion,\n(7) B* drag,         (8) Epoch,            (9) Satellite ID\n(0) Sat ID\nWhich parameter to change: ");
      scanf("%i",&i);
      if (i>=0 && i<=9) {
	printf("\nNew value: ");
	fgets(string,64,stdin);
	scanf("%s",string);
	if (i==0) strcpy(orb.desig,string);
	if (i==1) orb.eqinc=RAD(atof(string));
	if (i==2) orb.ascn=RAD(atof(string));
	if (i==3) orb.ecc=atof(string);
	if (i==4) orb.argp=RAD(atof(string));
	if (i==5) orb.mnan=RAD(atof(string));
	if (i==6) orb.rev=atof(string);
	if (i==7) orb.bstar=atof(string);
	if (i==8) {
	  orb.ep_year=2000+(int) floor(atof(string)/1000.0);
	  orb.ep_day=atof(string)-1000*floor(atof(string)/1000.0);
	}
	if (i==9) orb.satno=atoi(string);
	redraw=1;
	continue;
      }
      printf("\n================================================================================\n");
    }

    // Zoom
    if (c=='z') {
      click=1;
      mode=2;
    }

    // Execute zoom, or box delete
    if (c=='A') {
      if (click==0) {
	click=1;
      } else if (click==1 && mode==2) {
	xmin=(x0<x) ? x0 : x;
	xmax=(x0>x) ? x0 : x;
	ymin=(y0<y) ? y0 : y;
	ymax=(y0>y) ? y0 : y;

	click=0;
	mode=0;
	redraw=1;
	continue;
      } else {
	click=0;
	mode=0;
	redraw=1;
	continue;
      }
    }

    // Unzoom
    if (c=='r') {
      xmin=0.0;
      xmax=360.0;
      ymin=-90.0;
      ymax=90.0;
      mode=0;
      click=0;
      redraw=1;
      continue;
    }

    // Default tle
    if (c=='t') {
      orb.satno=d.p[0].satno;
      strcpy(orb.desig,d.p[0].desig);
      orb.ep_day=mjd2doy(0.5*(mjdmin+mjdmax),&orb.ep_year);
      satno=orb.satno;
      if (elset==0) {
	orb.eqinc=0.5*M_PI;
	orb.ascn=0.0;
	orb.ecc=0.0;
	orb.argp=0.0;
	orb.mnan=0.0;
	orb.rev=14.0;
	orb.bstar=0.5e-4;
	printf("LEO orbit\n");
      } else if (elset==1) {
	orb.eqinc=20.0*D2R;
	orb.ascn=0.0;
	orb.ecc=0.7;
	orb.argp=0.0;
	orb.mnan=0.0;
	orb.rev=2.25;
	orb.bstar=0.0;
	printf("GTO orbit\n");
      } else if (elset==2) {
	orb.eqinc=10.0*D2R;
	orb.ascn=0.0;
	orb.ecc=0.0;
	orb.argp=0.0;
	orb.mnan=0.0;
	orb.rev=1.0027;
	orb.bstar=0.0;
	printf("GSO orbit\n");
      } else if (elset==3) {
	orb.eqinc=63.4*D2R;
	orb.ascn=0.0;
	orb.ecc=0.7;
	orb.argp=0.0;
	orb.mnan=0.0;
	orb.rev=2.0;
	orb.bstar=0.0;
	printf("HEO orbit\n");
      }
      elset++;
      if (elset>3)
	elset=0;
      print_orb(&orb);
      printf("\n================================================================================\n");
      click=0;
      redraw=1;
      continue;
    }

    // Save
    x0=x;
    y0=y;
  }

  cpgend();

  free(d.p);

  fclose(stderr);
  return 0;
}

// Get observing site
struct site get_site(int site_id)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64];
  struct site s;
  char *env,filename[LIM];

  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/sites.txt",env);

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
    
    // Copy site
    if (id==site_id) {
      s.lat=lat;
      s.lng=lng;
      s.alt=alt;
      s.id=id;
      strcpy(s.observer,observer);
    }

  }
  fclose(file);

  return s;
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

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
}

// Observer position
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel)
{
  double ff,gc,gs,theta,s,dtheta;

  s=sin(lat*D2R);
  ff=sqrt(1.0-FLAT*(2.0-FLAT)*s*s);
  gc=1.0/ff+alt/XKMPER;
  gs=(1.0-FLAT)*(1.0-FLAT)/ff+alt/XKMPER;

  theta=gmst(mjd)+lng;
  dtheta=dgmst(mjd)*D2R/86400;

  pos->x=gc*cos(lat*D2R)*cos(theta*D2R)*XKMPER;
  pos->y=gc*cos(lat*D2R)*sin(theta*D2R)*XKMPER; 
  pos->z=gs*sin(lat*D2R)*XKMPER;
  vel->x=-gc*cos(lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
  vel->y=gc*cos(lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
  vel->z=0.0;

  return;
}

// Precess a celestial position
void precess(double mjd0,double ra0,double de0,double mjd,double *ra,double *de)
{
  double t0,t;
  double zeta,z,theta;
  double a,b,c;

  // Angles in radians
  ra0*=D2R;
  de0*=D2R;

  // Time in centuries
  t0=(mjd0-51544.5)/36525.0;
  t=(mjd-mjd0)/36525.0;

  // Precession angles
  zeta=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  zeta+=(0.30188-0.000344*t0)*t*t+0.017998*t*t*t;
  zeta*=D2R/3600.0;
  z=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  z+=(1.09468+0.000066*t0)*t*t+0.018203*t*t*t;
  z*=D2R/3600.0;
  theta=(2004.3109-0.85330*t0-0.000217*t0*t0)*t;
  theta+=-(0.42665+0.000217*t0)*t*t-0.041833*t*t*t;
  theta*=D2R/3600.0;
  
  a=cos(de0)*sin(ra0+zeta);
  b=cos(theta)*cos(de0)*cos(ra0+zeta)-sin(theta)*sin(de0);
  c=sin(theta)*cos(de0)*cos(ra0+zeta)+cos(theta)*sin(de0);

  *ra=(atan2(a,b)+z)*R2D;
  *de=asin(c)*R2D;

  if (*ra<360.0)
    *ra+=360.0;
  if (*ra>360.0)
    *ra-=360.0;

  return;
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
  if (year==1852 && month==10 && day<=4) b=0;

  jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

  return jd-2400000.5;
}

// Decode IOD Observations
struct point decode_iod_observation(char *iod_line)
{
  int year,month,iday,hour,min;
  int format,epoch,me,xe,sign;
  int site_id;
  double sec,ra,mm,ss,de,dd,ds,day,mjd0;
  char secbuf[6],sn[2],degbuf[3],buf1[3],buf2[6];
  struct point p;
  struct site s;
  xyz_t vel;

  // Strip newline
  iod_line[strlen(iod_line)-1]='\0';

  // Copy full line
  strcpy(p.iod_line,iod_line);

  // Set flag
  p.flag=1;

  // Get SSN
  sscanf(iod_line,"%5d",&p.satno);

  // Get desig
  sscanf(iod_line+6,"%s %s",buf1,buf2);
  sprintf(p.desig,"%s%s",buf1,buf2);

  // Get site
  sscanf(iod_line+16,"%4d",&site_id);
  s=get_site(site_id);

  // Decode date/time
  sscanf(iod_line+23,"%4d%2d%2d%2d%2d%5s",&year,&month,&iday,&hour,&min,secbuf);
  sec=atof(secbuf);
  sec/=pow(10,strlen(secbuf)-2);
  day=(double) iday+(double) hour/24.0+(double) min/1440.0+(double) sec/86400.0;
  p.mjd=date2mjd(year,month,day);

  // Get uncertainty in time
  sscanf(iod_line+41,"%1d%1d",&me,&xe);
  p.st=(float) me*pow(10,xe-8);

  // Get observer position
  obspos_xyz(p.mjd,s.lng,s.lat,s.alt,&p.obspos,&vel);

  // Skip empty observations
  if (strlen(iod_line)<64 || (iod_line[54]!='+' && iod_line[54]!='-')) 
    p.flag=0;

  // Get format, epoch
  sscanf(iod_line+44,"%1d%1d",&format,&epoch);

  // Read position
  sscanf(iod_line+47,"%2lf%2lf%3lf%1s",&ra,&mm,&ss,sn);
  sscanf(iod_line+55,"%2lf%2lf%2s",&de,&dd,degbuf);
  ds=atof(degbuf);
  if (strlen(degbuf)==1)
    ds*=10;
  sign=(sn[0]=='-') ? -1 : 1;
  sscanf(iod_line+62,"%1d%1d",&me,&xe);
  p.sr=(float) me*pow(10,xe-8);
  
  // Decode position
  switch(format) 
    {
      // Format 1: RA/DEC = HHMMSSs+DDMMSS MX   (MX in seconds of arc)
    case 1 : 
      ra+=mm/60+ss/36000;
      de=sign*(de+dd/60+ds/3600);
      p.sr/=3600.0;
      break;
      // Format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
    case 2:
      ra+=mm/60+ss/60000;
      de=sign*(de+dd/60+ds/6000);
      p.sr/=60.0;
      break;
      // Format 3: RA/DEC = HHMMmmm+DDdddd MX   (MX in degrees of arc)
    case 3 :
      ra+=mm/60+ss/60000;
      de=sign*(de+dd/100+ds/10000);
      break;
      // Format 7: RA/DEC = HHMMSSs+DDdddd MX   (MX in degrees of arc)
    case 7 :
      ra+=mm/60+ss/36000;
      de=sign*(de+dd/100+ds/10000);
      break;
    default :
      printf("IOD Format not implemented\n");
      p.flag=0;
      break;
  }
  // Convert to degrees
  ra*=15.0;

  // Get precession epoch
  if (epoch==0) {
    p.ra=ra;
    p.de=de;
    return p;
  } else if (epoch==4) {
    mjd0=33281.9235;
  } else if (epoch==5) {
    mjd0=51544.5;
  } else {
    printf("Observing epoch not implemented\n");
    p.flag=0;
  }

  // Precess position
  precess(mjd0,ra,de,p.mjd,&p.ra,&p.de);

  return p;
}

// Get a x and y from an AZI, ALT
void forward(double ra0,double de0,double ra,double de,double *x,double *y)
{
  int i;
  double phi,theta;
  struct celprm cel;
  struct prjprm prj;

  // Initialize Projection Parameters
  prj.flag=0;
  prj.r0=0.;
  for (i=0;i<10;prj.p[i++]=0.);

  // Initialize Reference Angles
  cel.ref[0]=ra0;
  cel.ref[1]=de0;
  cel.ref[2]=999.;
  cel.ref[3]=999.;
  cel.flag=0.;

  if (celset("STG",&cel,&prj)) {
    printf("Error in Projection (celset)\n");
    return;
  } else {
    if (celfwd("STG",ra,de,&cel,&phi,&theta,&prj,x,y)) {
      printf("Error in Projection (celfwd)\n");
      return;
    }
  }

  return;
}

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\t')
    c=' ';
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
}

// Read data
struct data read_data(char *filename)
{
  int i=0;
  char line[LIM];
  FILE *file;
  struct data d;

  // Open file
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"Failed to open %s\n",filename);
    exit(1);
  }

  // Count lines
  while (fgetline(file,line,LIM)>0) 
    i++;
  d.n=i;

  // Allocate
  d.p=(struct point *) malloc(sizeof(struct point)*d.n);

  // Rewind file
  rewind(file);

  // Read data
  i=0;
  while (fgetline(file,line,LIM)>0) 
    d.p[i++]=decode_iod_observation(line);

  // Close file
  fclose(file);

  return d;
}

// Chi-squared
double chisq(double a[])
{
  int i,imode,nsel;
  double chisq,rms;
  xyz_t satpos,satvel;
  double dx,dy,dz;
  double r;

  // Construct struct
  // a[0]: inclination
  // a[1]: RA of ascending node
  // a[2]: eccentricity
  // a[3]: argument of periastron
  // a[4]: mean anomaly
  // a[5]: revs per day

  if (a[2]<0.0)
    a[2]=0.0;
  if (a[0]<0.0) {
    a[0]*=-1;
    a[1]+=180.0;
  } else if (a[0]>180.0) {
    a[0]=180.0;
  }
  if (a[5]>20.0)
    a[5]=20.0;
  if (a[5]<0.1)
    a[5]=0.1;
  

  // Set parameters
  orb.eqinc=RAD(a[0]);
  orb.ascn=RAD(modulo(a[1],360.0));
  orb.ecc=a[2];
  orb.argp=RAD(modulo(a[3],360.0));
  orb.mnan=RAD(modulo(a[4],360.0));
  orb.rev=a[5];
  orb.bstar=a[6];

  // Initialize
  imode=init_sgdp4(&orb);
  if (imode==SGDP4_ERROR)
    printf("Error\n");

  // Loop over points
  for (i=0,nsel=0,chisq=0.0,rms=0.0;i<d.n;i++) {
    // Get satellite position
    satpos_xyz(d.p[i].mjd+2400000.5,&satpos,&satvel);
      
    // compute difference vector
    dx=satpos.x-d.p[i].obspos.x;  
    dy=satpos.y-d.p[i].obspos.y;
    dz=satpos.z-d.p[i].obspos.z;
      
    // Celestial position
    r=sqrt(dx*dx+dy*dy+dz*dz);
    d.p[i].rac=modulo(atan2(dy,dx)*R2D,360.0);
    d.p[i].dec=asin(dz/r)*R2D;
      
    // Compute offset
    forward(d.p[i].ra,d.p[i].de,d.p[i].rac,d.p[i].dec,&d.p[i].dx,&d.p[i].dy);
    d.p[i].dr=sqrt(d.p[i].dx*d.p[i].dx+d.p[i].dy*d.p[i].dy);

    if (d.p[i].flag==2) {
      // Compute chi-squared
      chisq+=pow(d.p[i].dr/d.p[i].sr,2);

      // Compute rms
      rms+=pow(d.p[i].dr,2);

      // Count selected points
      nsel++;
    }
  }
  if (nsel>0) 
    rms=sqrt(rms/(float) nsel);

  d.chisq=chisq;
  d.rms=rms;
  d.nsel=nsel;

  return chisq;
}

// Read tle
orbit_t read_tle(char *filename,int satno)
{
  int i;
  FILE *file;
  orbit_t orb;

  file=fopen(filename,"r");
  if (file==NULL) 
    fatal_error("Failed to open %s\n",filename);

  // Read TLE
  read_twoline(file,satno,&orb);
  fclose(file);

  return orb;
}

// MJD to DOY
double mjd2doy(double mjd,int *yr)
{
  int year,month,k=2;
  double day,doy;
  
  mjd2date(mjd,&year,&month,&day);

  if (year%4==0 && year%400!=0)
    k=1;

  doy=floor(275.0*month/9.0)-k*floor((month+9.0)/12.0)+day-30;

  *yr=year;

  return doy;
}

// Compute Date from Julian Day
void mjd2date(double mjd,int *year,int *month,double *day)
{
  double f,jd;
  int z,alpha,a,b,c,d,e;

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

  *day=b-d-floor(30.6001*e)+f;
  if (e<14)
    *month=e-1;
  else
    *month=e-13;

  if (*month>2)
    *year=c-4716;
  else
    *year=c-4715;

  return;
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
  sprintf(line2,"2 %05d %8.4f %8.4f %07.0f %8.4f %8.4f %11.8f    0",orb.satno,DEG(orb.eqinc),DEG(orb.ascn),1E7*orb.ecc,DEG(orb.argp),DEG(orb.mnan),orb.rev);

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

// Highlight
void highlight(float x0,float y0,float x,float y,int flag)
{
  int i;
  float xmin,xmax,ymin,ymax;

  xmin=(x0<x) ? x0 : x;
  xmax=(x0>x) ? x0 : x;
  ymin=(y0<y) ? y0 : y;
  ymax=(y0>y) ? y0 : y;
  for (i=0;i<d.n;i++) 
    if (d.p[i].ra>xmin && d.p[i].ra<xmax && d.p[i].de>ymin && d.p[i].de<ymax && d.p[i].flag!=0)
      d.p[i].flag=flag;

  return;
}

// Select time range
void time_range(double *mjdmin,double *mjdmax,int flag)
{
  int i,n;
  float c;

  for (i=0,n=0;i<d.n;i++) {
    if (d.p[i].flag==flag) {
      if (n==0) {
	*mjdmin=d.p[i].mjd;
	*mjdmax=d.p[i].mjd;
      }
      if (d.p[i].mjd< *mjdmin) *mjdmin=d.p[i].mjd;
      if (d.p[i].mjd> *mjdmax) *mjdmax=d.p[i].mjd;
      n++;
    }
  }
  c=0.1*(*mjdmax- *mjdmin);
  *mjdmin-=c;
  *mjdmax+=c;

  return;
}

// Print TLE
void print_tle(orbit_t orb,char *filename)
{
  int i,n;
  FILE *file;
  double mjdmin,mjdmax;
  int year,month;
  double day;
  char line1[70],line2[70];

  // Count number of points
  for (i=0,n=0;i<d.n;i++) {
    if (d.p[i].flag==2) {
      if (n==0) {
	mjdmin=d.p[i].mjd;
	mjdmax=d.p[i].mjd;
      }
      if (d.p[i].mjd<mjdmin) mjdmin=d.p[i].mjd;
      if (d.p[i].mjd>mjdmax) mjdmax=d.p[i].mjd;
      n++;
    }
  }

  // Write TLE
  file=fopen(filename,"w");
  format_tle(orb,line1,line2);
  fprintf(file,"OBJ\n%s\n%s\n",line1,line2);

  mjd2date(mjdmin,&year,&month,&day);
  fprintf(file,"# %4d%02d%05.2lf-",year,month,day);
  mjd2date(mjdmax,&year,&month,&day);
  fprintf(file,"%4d%02d%05.2lf, %d measurements, %.3lf deg rms\n",year,month,day,n,d.rms);
  fclose(file);

  return;
}

// Fit
void fit(orbit_t orb,int *ia)
{
  int i,n;
  double a[7],da[7];
  //  double db[7]={5.0,5.0,0.1,5.0,5.0,0.5,0.0001};
  double db[7]={1.0,1.0,0.02,1.0,1.0,0.1,0.0001};

  a[0]=orb.eqinc*R2D;
  da[0]=da[0]*R2D;
  a[1]=orb.ascn*R2D;
  da[1]=da[1]*R2D;
  a[2]=orb.ecc;
  a[3]=orb.argp*R2D;
  da[3]=da[3]*R2D;
  a[4]=orb.mnan*R2D;
  da[4]=da[4]*R2D;
  a[5]=orb.rev;
  a[6]=orb.bstar;

  for (i=0;i<7;i++) {
    if (ia[i]==1)
      da[i]=db[i];
    else
      da[i]=0.0;
  }

  // Construct struct
  // a[0]: inclination
  // a[1]: RA of ascending node
  // a[2]: eccentricity
  // a[3]: argument of periastron
  // a[4]: mean anomaly
  // a[5]: revs per day
  // a[6]: bstar

  // Count highlighted points
  for (i=0,n=0;i<d.n;i++)
    if (d.p[i].flag==2)
      n++;

  if (n>0)
    versafit(n,7,a,da,chisq,0.0,1e-6,"n");

  // Return parameters
  orb.eqinc=RAD(a[0]);
  orb.ascn=RAD(modulo(a[1],360.0));
  orb.ecc=a[2];
  orb.argp=RAD(modulo(a[3],360.0));
  orb.mnan=RAD(modulo(a[4],360.0));
  orb.rev=a[5];
  orb.bstar=a[6];

  return;
}

void usage()
{
  printf("satfit d:c:i:haCo:p\n\n");
  printf("d    IOD observations\n");
  printf("c    TLE catalog\n");
  printf("i    Satellite ID (NORAD)\n");
  printf("C    Fit circular orbit\n");
  printf("p    No plotting\n");
  printf("o    Output TLE file\n");
  printf("a    Adjust MA and RAAN\n");
  printf("h    This help\n");

  return;
}
