/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
#include "grid/cartesian1D.h"
#include "./saint-venantNN.h"

// Rheological properties
#define tauY 3.4765171
#define muN 0.10
#define nHB 1.00

#define rhoFluid 1000.0
#define gravity 9.81
#define sinTheta 0.06
#define cosTheta (pow(1.0-sinTheta*sinTheta,0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta)

#define froude 0.97601
#define normalDepth 0.019688057
#define normalVel (froude*pow(gPrime*normalDepth,0.50))
#define normalMaxVel (normalVel/(1.0-nHB/(2.0*nHB+1.0)*(1.0-alphaCoeff)))
#define colorbarMax (normalMaxVel*2.50)

#define alphaCoeff (tauY/(rhoFluid*gravity*normalDepth*sinTheta))

#define piVal 3.14159265
#define disAmp 0.05

#define DOMAINLENGTH (2.0*piVal/1.2*normalDepth/So)
#define PLOTRANGEMAX (2.01*normalDepth)

#define INITLEVEL 10

#define simTime 35.0
#define outputInterval 0.50

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

// double nu_eq(double shear,double pipi){
double nu_eq(double shear){
  double nu_eq;
  nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)));
  // nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)))*(1.-exp((-1.)*papaCoeff*sqrt(sq(shear) + sq(1.e-10))));
  return nu_eq;
}

double hbProfile (double z, double ht)
{
  // general normal-flow profile of Herschel-Bulkley fluids
  double zo;
  zo = ht*(1.-alphaCoeff);
  // zo = normalDepth*alphaCoeff;
  // if (z<=(normalDepth*alphaCoeff)) {
  if (z<=zo) {
    // shear zone
    return (nHB/(nHB+1.)*pow((rhoFluid*gravity*sinTheta/muN*pow(zo, (nHB+1.))), (1./nHB))*(1.-pow((1.-z/zo),((nHB+1.)/nHB))));
  }
  else {
    // plug zone
    return (nHB/(nHB+1.)*pow((rhoFluid*gravity*sinTheta/muN*pow(zo, (nHB+1.))), (1./nHB)));
  }
}

char s[80];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = gPrime;
  init_grid(1 << (INITLEVEL));
  nl = 32;
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
    double totalDepth = normalDepth*(1.0+disAmp*sin(2. * pi * x / DOMAINLENGTH));
    double z = zb[];
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = totalDepth;
    for (int l = 0; l < nl; l++) {
      z += h[]/2.;
      u = ul[l];
      u.x[] = hbProfile(z, totalDepth);
      z += h[]/2.;
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
    boundary ((scalar *){u});
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval){
  sprintf (s, "slice-%g.txt", t);
  FILE * fp1 = fopen (s, "w"); 
  foreach(serial){
    double zCoord = zb[];
    for (int l = 0; l < nl; l++) {
      zCoord += layer[l]*h[]*0.50;
      u = ul[l];
      fprintf (fp1, "%g %g %g \n", x, zCoord, u.x[] );
      zCoord += layer[l]*h[]*0.50;
    }
  }
  fclose(fp1);

  sprintf (s, "depth-%g.txt", t);
  FILE * fp2 = fopen (s, "w"); 
  foreach (serial) {
    fprintf (fp2, "%g %g \n", x, h[]);
  }
  fclose(fp2);
}

/**
Post-processing visualization module
*/ 
void setup (FILE * fp)
{
  // FIXME: other customized color schemes
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [0:%g]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [%g:%g]\n"
	   "set yrange [0.0:%g]\n"
     "set lmargin at screen 0.1\n"
     "set rmargin at screen 0.9\n", colorbarMax, 0.0, DOMAINLENGTH, PLOTRANGEMAX
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %g'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    // fprintf (fp, "%g %g %g\n", x, max(0.5, z), u.x[]);
    for (int l = 0; l < nl; l++) {
      z += layer[l]*h[]*0.50;
      u = ul[l];
      // z += h[]/2.0;
      fprintf (fp, "%g %g %g\n", x, z, u.x[]);
      // z += h[]/2.0;
      z += layer[l]*h[]*0.50;
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += outputInterval; t <= simTime)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,320\n"
	   "set output 'plot-%g.png'\n", t);
  if (i == 0)
    setup (fp);
  plot (fp);
}

// event moviemaker (t = end) {
//     system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
// 	    "ppm2mp4 movie.mp4");
// }