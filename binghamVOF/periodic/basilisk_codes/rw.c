/**
# 2D roll wave simulation

All variables are in SI units
*/

#include "navier-stokes/centered.h"
#define FILTERED 1.0

#include "two-phaseVP.h"

// Domain extent
#define LDOMAIN 20.926
#define CHANNELSIN 0.0680335
#define CHANNELCOS pow(1.0-pow(CHANNELSIN,2.0),0.50)
#define GRAVCOEFF 9.81

// normal-flow definition
double  normalDepth = 0.237831;
double normalVelocitySurface = 0.0;

// disturbance parameters
double distAmp = 0.225;
// if use non-periodic BC, a wavelength is also needed

// Maximum refinement: guarantee more than 23 cell through normal depth
#define MAXLEVEL 11
#define MINLEVEL 8
#define INITLEVEL 8


// max run-time
double tmax = 29.0;

double tsnap = 4.0;

// mesh adaptivity parameters
double uemax = 0.01;


// steady uniform flow solution (parabolic velocity profile)
double velDist( double y, double x, double YS, double rhogx, double HN, double mu0, double Lx, double amp){
    double yieldY = HN-YS/(rhogx);
    double Up = pow((YS-HN*(rhogx)),2.0)/(2.0*mu0*(rhogx));
    double hDist = HN*(1.0+amp*sin(2.0*pi*x/Lx));
    return y < hDist ? (y<yieldY ? (Up*(1.0-pow(((y-yieldY)/yieldY),2.0))) : (Up)) : 0.0;
}

// initial VOF field setting
double heightDist( double x,  double HN, double Lx, double amp){
    return (HN*(1.0+amp*sin(2.0*pi*x/Lx)));
}

// hydrostatic pressure distribution
double hydroPressureDist( double y,  double HN, double gy, double rhoFluid){
    return (y<HN ? (rhoFluid*gy*y) : 0.0);
}

/**
Boundary conditions */

// p[top] = dirichlet(-RHOF*LDOMAIN);
// pf[top] = dirichlet(-RHOF*LDOMAIN);
p[top] = dirichlet(0.0);
pf[top] = dirichlet(0.0);
u.n[top] = neumann(0);
u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);

/**
the three cases in the main */

int main() 
{
  size (LDOMAIN);
  // number of grid points
  // initial grid is coarse to save memory. to be refined locally
  N = 1 << INITLEVEL;
  // maximum timestep
  DT = 1e-3;
  TOLERANCE = 1e-3;

  // assign material properties
  // 1: fluid, 2: gas
  rho1 = 2.10e3; mu1 = 10.0; rho2 = 1.12; mu2 = 1e-6;

  tauy = 100.0;
  mumax = 2.0e4;

  periodic (right);

  const face vector g[] = {CHANNELSIN*GRAVCOEFF,(-1.0)*GRAVCOEFF*CHANNELCOS};
  a = g;
  // by default, the origin is the lower-left corner
  // origin(...,...)


  run();
}

/** initial conditions
*/
event init (i = 0) 
{
  refine((y < normalDepth*2.25) && level < MAXLEVEL);
  CFL = 0.3;

  scalar phi[];
  foreach_vertex()
    phi[] = heightDist(x, normalDepth, LDOMAIN, distAmp)-y;
  // pay attention to the definition of heavy phase and light phase!
  fractions (phi, f);
/*
initialisation of hydrostatic pressure for heacy phase  
*/

  foreach() {
    u.x[] = velDist( y, x, tauy, (GRAVCOEFF*rho1*CHANNELSIN), normalDepth, mu1, LDOMAIN, distAmp);
    u.y[] = 0.0;
    p[]= hydroPressureDist(y, normalDepth, (CHANNELCOS*GRAVCOEFF), (rho1));
  }
}

// event timingParameters (i = 0; i+=20) 
// {
//   static FILE * fp = fopen ("parameters", "w");
//   fprintf (fp, "T = %g\tDT = %g\tCFL = %g\n", t, dt, CFL);
//   fclose(fp); 
// }

// event timingsOutput (i+=10) 
// {
//   fprintf(stdout,"%g\t%g\t%d\t\n",t,dt,i);
// }

event timingLogfile (i += 10)
{
  fprintf (stderr, "%g %g \n", t, dt);
}

/**
txt files outputs
*/
event txtFileFields (t = 0.0; t += 0.2; t <= tmax) 
{
  char s[50], s1[50];
  sprintf (s, "freeSurface-%g.txt", t);
  FILE * fpf = fopen (s, "w");
  output_facets (f, fpf); 
  fclose (fpf);

  sprintf (s1, "field-%g.txt", t);
  FILE * fptxt = fopen (s1, "w");
  output_field ({f,p,u,uf,pf}, fptxt, linear = true);
  fclose (fptxt);
}

/**
gfs files output
*/
event gfsSnapshots (t = 0.0; t += 0.1; t <= tmax)
{
  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  output_gfs (file = name, t = t, list = {u.x,u.y,f});
}

event writingFiles (t += tsnap; t <= tmax) {
  char nameOut[80];
  dump (file = "dump");
  sprintf (nameOut, "dumpSnapshot-%6.5f", t);
  dump(file=nameOut);
}

/**
 * TODO:Adaptivity
*/
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.0025,uemax,uemax,uemax}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

// event adapt (i++) {
//   adapt_wavelet ({f}, (double[]){0.01}, MAXLEVEL, MINLEVEL);
// }
