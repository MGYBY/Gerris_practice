/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
#include "grid/cartesian1D.h"
#include "./saint-venant-cb.h"

// Rheological properties
#define tauY 3.4765171
#define muN 0.10

#define rhoFluid 1000.0
#define nuVal (muN/rhoFluid)
#define gravity 9.81
#define sinTheta 0.06
#define cosTheta (pow(1.0-sinTheta*sinTheta,0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta)

#define froude 0.97601
#define normalDepth 0.019688057
#define normalVel (froude*pow(gPrime*normalDepth,0.50))

#define xiNormParam (tauY/(rhoFluid*gravity*normalDepth*sinTheta))
#define alpha1Param(xiParam) ((1.-xiParam)*(1.+xiParam/2.))
#define alpha2Param(xiParam) ((1.-xiParam)*(1.+5./4.*xiParam)*pow((1.+xiParam/2.),(-2.0)))
#define beta1Param(xiParam) (xiParam*pow((1.-xiParam),5.0)*(1.+5./4.*xiParam)*(((1.+5./2.*xiParam+pow(xiParam,2.)/2.)*(1.+2.*xiParam+3.*pow(xiParam,2.)+3./2.*pow(xiParam,3.0)))/(1.+61./16.*xiParam+63./8.*pow(xiParam,2.)+63./8.*pow(xiParam,3.)+21./16.*pow(xiParam,4.))))
#define beta2Param(xiParam) ((1.-xiParam)*(1.+5./4.*xiParam)*((1.+2.*xiParam+3.*pow(xiParam,2.0)+3./2.*pow(xiParam,3.0))/(1.+61./16.*xiParam+63./8.*pow(xiParam,2.0)+63./8.*pow(xiParam,3.)+21./16.*pow(xiParam,4.))))
#define rParam(xiParam) ((1.+11./4.*xiParam+87./16.*pow(xiParam,2.)+107./16.*pow(xiParam,3.)+pow(xiParam,4.))/((1.+xiParam/2.)*(1.+2.*xiParam+3.*pow(xiParam,2.)+3./2.*pow(xiParam,3.))))

#define normalPhi ((pow(normalVel,2.0)/(5.0*pow(normalDepth,2.0)))*alpha2Param(xiNormParam))

#define piVal 3.14159265
#define disAmp 0.10

#define DOMAINLENGTH (2.0*piVal/1.2*normalDepth/So)

#define INITLEVEL 10

#define simTime 30.0
#define outputInterval 0.50

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

char s[80];
scalar depthGrad[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = gPrime;
  init_grid(1 << (INITLEVEL));

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
    zb[] = 0.0;
    h[] = totalDepth;
    // simple velocity field for now
    u.x[] = normalVel;
    hie[] = 0.5*(pow(normalVel,2.0)+pow(normalDepth,2.0)*normalPhi+gravity*cosTheta*normalDepth)*normalDepth;

    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    }
}

static double momFric(double u, double h, double hie)
{
     double rhs, xiVal, phiVal;
     xiVal = tauY/(rhoFluid*gravity*h*sinTheta);
     phiVal = h>dry ? (2.*hie-h*(gPrime*h+pow(u,2.0)))/pow(h,3.0) : 0.0;
     rhs = h>dry ? (gPrime*So-tauY/(rhoFluid*h)-3.*nuVal*u/(pow(h,2.0)*alpha1Param(xiVal)))*(alpha1Param(xiVal)+7./360.*((pow((gPrime*h*So),2.0))/(pow(nuVal,2.0)*phiVal))*beta1Param(xiVal)) + 7./6.*(gPrime*So/phiVal)*(phiVal-pow(u,2.0)/(5.*pow(h,2.0))*alpha2Param(xiVal))*beta2Param(xiVal) : 0.0;
     return rhs;
}

static double hieFric(double u, double h, double hie)
{
     double rhs, xiVal, phiVal;
     xiVal = tauY/(rhoFluid*gravity*h*sinTheta);
     phiVal = h>dry ? (2.*hie-h*(gPrime*h+pow(u,2.0)))/pow(h,3.0) : 0.0;
     rhs = h>dry ? u*h*(gPrime*So-tauY/(rhoFluid*h)-3.*nuVal*u/(pow(h,2.0)*alpha1Param(xiVal)))*(alpha1Param(xiVal)+7./1080.*((pow((gPrime*h*So),2.0))/(pow(nuVal,2.0)*phiVal))*beta1Param(xiVal)*rParam(xiVal)) + 7./18.*(u*gPrime*h*So/phiVal)*(phiVal-pow(u,2.0)*alpha2Param(xiVal)/(5.*pow(h,2.0)))*beta2Param(xiVal)*rParam(xiVal) : 0.0;
     return rhs;
}

/**
  Source term.
*/
event frictionTerm (i++) {
  foreach(){
    double uPrev, uMed, hieMed;
    uPrev = u.x[];

    // x-mom eqn source term
    uMed = u.x[] + dt *  momFric(u.x[], h[], hie[]);
    uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * momFric(uMed, h[], hie[]);
    u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * momFric(uMed, h[], hie[]);

    // internal energy eqn he source term
    hieMed = hie[] + dt *  hieFric(uPrev, h[], hie[]);
    hieMed = (3. / 4.) * hie[]+ (1. / 4.) * hieMed + (1. / 4.) * dt * hieFric(uPrev, h[], hieMed);
    hie[] = (1. / 3.) * hie[] + (2. / 3.) * hieMed + (2. / 3.) * dt * hieFric(uPrev, h[], hieMed);

    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    }

    boundary ((scalar *){u, hie});
}

/*
 * AMR here
 *
 */
// event adapt1 (i++) {
//   adapt_wavelet({h, depthGrad, u.x[]}, (double[]){normalDepth/500.0, 0.004, normalVel/400.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
// //      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
// //   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
// }

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval) {
  sprintf (s, "slice-%g.txt", t);
  FILE * fp2 = fopen (s, "w"); 
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g \n", x, h[], u.x[], hie[]);
  }
  fclose(fp2);
}

/**
Visualise the wave profile
*/

void plot_profile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't%g.png'\n"
             "set title 't = %.2f'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
               fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplotOutput1(t = 0; t <= simTime; t += outputInterval)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
}

// event moviemaker (t = end) {
//     system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
// 	    "ppm2mp4 movie.mp4");
// }
