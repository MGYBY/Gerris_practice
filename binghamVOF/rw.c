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
#define GRAVCOEFF 9.81

// normal-flow definition
double  normalDepth = 0.237831;
double normalVelocitySurface = 0.0;
double channelCos = pow(1.0-pow(CHANNELSIN,2.0),0.50);

// disturbance parameters
double distAmp = 0.0;
// if use non-periodic BC, a wavelength is also needed

// Maximum refinement: guarantee more than 23 cell through normal depth
#define MAXLEVEL 11
#define MINLEVEL 8
#define INITLEVEL 4


// max run-time
double tmax = 0.059;

// mesh adaptivity parameters
// double uemax = ;

// a file containing the free-surface
FILE * fpf;


// steady uniform flow solution (parabolic velocity profile)
double velDist( double y, double YS, double rhogx, double HN, double mu0){
    double yieldY = HN-YS/(rhogx);
    double Up = pow((YS-HN*(rhogx)),2.0)/(2.0*mu0*(rhogx));
    return (y<yieldY ? (Up*(1.0-pow(((y-yieldY)/yieldY),2.0))) : (Up));
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

int main() {
  L0 = LDOMAIN;
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

  const face vector g[] = {CHANNELSIN*GRAVCOEFF,(-1.0)*GRAVCOEFF*channelCos};
  a = g;
  // by default, the origin is the lower-left corner
  // origin(...,...)

  fpf = fopen ("freeSurface", "w");
  fclose (fpf);

  run();
}
/**
initial heap, a rectangle
*/

event init (t = 0) {

  refine((y<=normalDepth) && level < MAXLEVEL);

  scalar phi[];
  foreach_vertex()
    phi[] = heightDist(x, normalDepth, LDOMAIN, distAmp)-y;
  // pay attention to the definition of heavy phase and light phase!
  fractions (phi, f);
/*
initialisation of hydrostatic pressure for heacy phase  
*/

  foreach() {
    u.x[] = velDist(y, tauy, (rho1*GRAVCOEFF*CHANNELSIN), normalDepth, mu1);
    u.y[] = 0.0;
    p[]= hydroPressureDist(y, normalDepth, (channelCos*GRAVCOEFF), (rho1));
  }
}

event timingsOutput (i+=10) {
  fprintf(stdout,"%g\t%g\t%d\t\n",t,dt,i);
}

/**
txt files outputs
*/
// event txtFileFields (t = 0.0; t += 0.05; t <= tmax) {
//   output_facets (f, fpf); 
//   char s[80];
//   sprintf (s, "field-%g.txt", t);
//   FILE * fptxt = fopen (s, "w");
//   output_field ({f,p,u,uf,pf}, fptxt, linear = true);
//   fclose (fptxt);
// }

/**
gfs files output
*/
event gfsSnapshots (t = 0.0; t += 0.001; t <= tmax) {
  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  output_gfs (file = name, t = t, list = {u.x,u.y,f});
}

/**
 * TODO:Adaptivity
*/
// event adapt (i++) {
//   adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
// }

// event adapt (i++) {
//   adapt_wavelet ({f}, (double[]){0.01}, MAXLEVEL, MINLEVEL);
// }