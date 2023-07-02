/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
#include "grid/bitree.h"
// #include "grid/cartesian1D.h"
#include "./saint-venantNN.h"

// Rheological properties
#define tauY 0.0
#define muN 150.0
#define nHB 0.30

#define rhoFluid 2130.0
#define gravity 9.81
#define sinTheta 0.00
#define cosTheta (pow(1.0-sinTheta*sinTheta,0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta)

#define aspectRatio 1.00
#define initDepth 1.00

#define colorbarMax (pow((gPrime*initDepth), 0.50)*0.50)

#define simTime 3000.0

// #define alphaCoeff (tauY/(rhoFluid*gravity*normalDepth*sinTheta))

#define piVal 3.14159265
#define disAmp 0.10

// initial hump
#define humpXCenter (2.0*normalDepth/So)
#define humpWidth (1.0*normalDepth/So)

#define OUTPUTTIME1 0.2
#define OUTPUTTIME2 2.00
#define OUTPUTTIME3 50.00
#define OUTPUTPERIOD 20.0
#define OUTPUTPERIOD2 200.0
#define OUTPUTPLOT 10.00
#define DOMAINLENGTH (initDepth*40.0)
#define PLOTRANGEMAX (1.05*initDepth)

#define INITLEVEL 13
#define MAXLEVEL 13
#define MINLEVEL 4

#define XFROPERIOD 4.0

#define initialStageTime 3.0
#define restrictDt 1.0e-4

#define inletRefLen (aspectRatio*3.0)

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

// double nu_eq(double shear,double pipi){
double nu_eq(double shear){
  double nu_eq;
  nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.1e-8)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.1e-8)));
  // nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)))*(1.-exp((-1.)*papaCoeff*sqrt(sq(shear) + sq(1.e-10))));
  return nu_eq;
}

// double hbProfile (double z, double ht)
// {
//   // general normal-flow profile of Herschel-Bulkley fluids
//   double zo;
//   zo = ht*(1.-alphaCoeff);
//   // zo = normalDepth*alphaCoeff;
//   // if (z<=(normalDepth*alphaCoeff)) {
//   if (z<=zo) {
//     // shear zone
//     return (nHB/(nHB+1.)*pow((rhoFluid*gravity*sinTheta/muN*pow(zo, (nHB+1.))), (1./nHB))*(1.-pow((1.-z/zo),((nHB+1.)/nHB))));
//   }
//   else {
//     // plug zone
//     return (nHB/(nHB+1.)*pow((rhoFluid*gravity*sinTheta/muN*pow(zo, (nHB+1.))), (1./nHB)));
//   }
// }

char s[80];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = gPrime;
  init_grid(1 << (INITLEVEL));
  nl = 50;
  nu = 1.; // dummy

  CFL = 0.40;

//   periodic (right);

  run();
}

scalar depthGrad[];
scalar fr[], frGrad[], velAve[];

/**
## Initialization  */
event init (i = 0) {
  /**
  We initialize *h*. */
  foreach(){
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = x<=(initDepth*aspectRatio) ? initDepth : 0.0;
    depthGrad[] = fabs(h[-1]-h[1])/Delta;
    velAve[] = 0.0;

    for (int l = 0; l < nl; l++) {
      u = ul[l];
      u.x[] = 0.0;
      velAve[] += h[]>dry ? u.x[]*h[]*layer[l] : 0.0;
    }
    velAve[] = h[]>dry ? velAve[]/h[] : 0.0;
    fr[] = h[]>dry ? velAve[]/pow((G*h[]), 0.50) : 0.0;
  }
}

event dtMax(t = 0; t<=initialStageTime; t+=restrictDt)
{

}

  /**
## Output
  We print the elevation  */
event acceleration (i++) {
  foreach(){
    depthGrad[] = fabs(h[-1]-h[1])/Delta;
    velAve[] = 0.0;

    for (int l = 0; l < nl; l++) {
      u = ul[l];
      u.x[] += dt*gravity*sinTheta;
      velAve[] += h[]>dry ? u.x[]*h[]*layer[l] : 0.0;
    }
    velAve[] = h[]>dry ? velAve[]/h[] : 0.0;
    fr[] = h[]>dry ? velAve[]/pow((G*h[]), 0.50) : 0.0;
    frGrad[] = fabs(fr[-1]-fr[1])/Delta;
    }
    boundary ((scalar *){u});
}

event hmaxUmax1 (t=0; i+=4; t<XFROPERIOD)
{
     double maxDepth = 0.0;
     double maxVel = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     double maxVelLocX = 0.0;
     double maxVelDepth = 0.0;
     FILE *fp2 = fopen("maxDepth", "a+");
     FILE *fp3 = fopen("maxVel", "a+");
     FILE *fp4 = fopen("maxFr", "a+");
     foreach ()
     {
           if (h[]>maxDepth){
                maxDepth = h[];
                maxDepthLocX = x;
           }

           if (velAve[]>maxVel && h[]>dry){
                maxVel = velAve[];
                maxVelLocX = x;
           }
     }
//      stats s1 = statsf (re);
     stats s2 = statsf (fr);
     fprintf(fp2, "%.10g %.10g %.10g \n", t, maxDepth, maxDepthLocX);
     fprintf(fp3, "%.10g %.10g %.10g \n", t, maxVel, maxVelLocX);
     fprintf(fp4, "%.10g %.10g \n", t, s2.max);
     fclose(fp2);
     fclose(fp3);
     fclose(fp4);
}

event hmaxUmax2 (t=XFROPERIOD; t<=simTime; i+=50)
{
     double maxDepth = 0.0;
     double maxVel = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     double maxVelLocX = 0.0;
     double maxVelDepth = 0.0;
     FILE *fp2 = fopen("maxDepth", "a+");
     FILE *fp3 = fopen("maxVel", "a+");
     FILE *fp4 = fopen("maxFr", "a+");
     foreach ()
     {
           if (h[]>maxDepth){
                maxDepth = h[];
                maxDepthLocX = x;
           }

           if (velAve[]>maxVel && h[]>dry){
                maxVel = velAve[];
                maxVelLocX = x;
           }
     }
//      stats s1 = statsf (re);
     stats s2 = statsf (fr);
     fprintf(fp2, "%.10g %.10g %.10g \n", t, maxDepth, maxDepthLocX);
     fprintf(fp3, "%.10g %.10g %.10g \n", t, maxVel, maxVelLocX);
     fprintf(fp4, "%.10g %.10g \n", t, s2.max);
     fclose(fp2);
     fclose(fp3);
     fclose(fp4);
}

event hFront1 (t=0; i+=4; t<XFROPERIOD)
{
     double aveDepth = 0.0;
     double aveVel = 0.0;
     double xf = 0.0;
     FILE *fp3 = fopen("frontPosMod", "a+");
     foreach ()
     {
//           frontXPos[] = h[]>=dry ? x : 0.0;
          xf = h[] > dry ?  max(xf,x) :  xf ;
          aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
          aveVel += (h[]>=dry ? Delta*velAve[] : 0.0);
     }
     stats s1 = statsf (h);
//      stats s2 = statsf (frontXPos);
//      xf = s2.max;
     fprintf(fp3, "%.10g %.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf), (s1.max/xf));
     fclose(fp3);
}

event hFront2 (t=XFROPERIOD; t<=simTime; i+=80)
{
     double aveDepth = 0.0;
     double aveVel = 0.0;
     double xf = 0.0;
     FILE *fp3 = fopen("frontPosMod", "a+");
     foreach ()
     {
//           frontXPos[] = h[]>=dry ? x : 0.0;
          xf = h[] > dry ?  max(xf,x) :  xf ;
          aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
          aveVel += (h[]>=dry ? Delta*velAve[] : 0.0);
     }
     stats s1 = statsf (h);
//      stats s2 = statsf (frontXPos);
//      xf = s2.max;
     fprintf(fp3, "%.10g %.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf), (s1.max/xf));
     fclose(fp3);
}

event adapt1 (i++) {
  astats s = adapt_wavelet({h, depthGrad, frGrad}, (double[]){initDepth/300.0, 2.e-5, 3.5e-4}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output1  (t = 0.0; t <= OUTPUTPERIOD; t += OUTPUTTIME1){
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
    fprintf (fp2, "%g %g %g \n", x, h[], velAve[]);
  }
  fclose(fp2);
}

event output2  (t = OUTPUTPERIOD; t <= OUTPUTPERIOD2; t += OUTPUTTIME2){
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
    fprintf (fp2, "%g %g %g \n", x, h[], velAve[]);
  }
  fclose(fp2);
}

event output3  (t = OUTPUTPERIOD2; t <= simTime; t += OUTPUTTIME3){
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
    fprintf (fp2, "%g %g %g \n", x, h[], velAve[]);
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

event gnuplot (t += OUTPUTPLOT; t <= simTime)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set term pngcairo font \",10\" size 683,213\n"
	   "set output 'plot-%g.png'\n", t);
  if (i == 0)
    setup (fp);
  plot (fp);
}

// event moviemaker (t = end) {
//     system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
// 	    "ppm2mp4 movie.mp4");
// }
