/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
#include "grid/cartesian1D.h"
#include "./saint-venantNN.h"

// Rheological properties
#define tauY 2.84201
#define muN 0.10
#define nHB 1.00

#define rhoFluid 1000.0
#define gravity 9.81
#define sinTheta 0.06
#define cosTheta (pow(1.0-sinTheta*sinTheta,0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta)

#define froude 0.7214
#define normalDepth 0.016094745
#define normalVel (froude*pow(gPrime*normalDepth,0.50))

#define DOMAINLENGTH (20.0*normalDepth)

#define simTime 50.0

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

// double nu_eq(double shear,double pipi){
double nu_eq(double shear){
  double nu_eq;
  nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)));
  return nu_eq;
}
/**



*/

char s[80];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = gPrime;
  N  = 8;
  nl = 75;
  nu = 1.; // dummy

  CFL = 0.40;

  periodic (right);

  run();
}

/**
## Initialization  */
event init (i = 0) {
  /**
  We initialize *h*. */
  foreach(){
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = normalDepth;
    for (int l = 0; l < nl; l++) {
      u = ul[l];
      u.x[] = 0.0;
    }
    }
}

  /**
## Output
  We print the elevation  */
event acceleration (i++) {
  foreach(){
    for (int l = 0; l < nl; l++) {
      u = ul[l];
      u.x[] += dt*gravity*sinTheta;
    }
    }
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=1.0){
  sprintf (s, "slice-%g.txt", t);
  FILE * fp = fopen (s, "w"); 
  foreach(){
    double zCoord = 0.;
    for (int l = 0; l < nl; l++) {
      zCoord += layer[l]*h[];
      u = ul[l];
      fprintf (fp, "%g %g %g \n", x, zCoord, u.x[] );
    }
  }
  fclose(fp);
}
