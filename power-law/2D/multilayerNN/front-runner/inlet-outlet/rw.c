/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
// #include "grid/cartesian1D.h"
#include "grid/bitree.h"
#include "./saint-venantNN.h"

// Rheological properties
#define bParam 0.30
#define nHB 0.30
#define froude 1.50
#define grav (1.0/pow(froude,2.0))
#define rParam (pow(((nHB/(nHB+1.0))*pow(grav,(1.0/nHB))*pow((1.0-bParam),((nHB+1.0)/nHB))*(1.0-(nHB/(2.0*nHB+1.0))*(1.0-bParam))), nHB))

#define normalDepth 1.00
// #define normalDepth 0.095263786
#define normalVel 1.00
// #define alphaCoeff (tauY/(rhoFluid*gravity*normalDepth*sinTheta))
#define normalMaxVel (normalVel/(1.0-nHB/(2.0*nHB+1.0)*(1.0-bParam)))
#define colorbarMax (normalMaxVel*3.50)

#define piVal 3.14159265
#define disAmp 0.10

#define DOMAINLENGTH (100.0)
#define PLOTRANGEMAX (5.00*normalDepth)

#define distPeriod 2.0
#define refLen (distPeriod*2.50)

#define MAXLEVEL 14
#define MINLEVEL 3
#define INITLEVEL 14

#define simTime 32.0
#define outputInterval 2.00

#define visRegThre 1.01e-10
#define yieldSurfThre 2.e-7

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

// double nu_eq(double shear,double pipi){
double nu_eq(double shear){
  double nu_eq;
  nu_eq = rParam*pow(sqrt(sq(shear) + sq(visRegThre)), (nHB-1.0)) + (bParam/pow(froude,2.0))*1.0/(sqrt(sq(shear) + sq(visRegThre)));
  // nu_eq = muN/rhoFluid*pow(sqrt(sq(shear) + sq(1.e-10)), (nHB-1.0))+tauY/rhoFluid/(sqrt(sq(shear) + sq(1.e-10)))*(1.-exp((-1.)*papaCoeff*sqrt(sq(shear) + sq(1.e-10))));
  return nu_eq;
}

double hbProfile (double z, double ht)
{
  // general normal-flow profile of Herschel-Bulkley fluids
  double zo;
  zo = ht*(1.-bParam);
  // zo = normalDepth*alphaCoeff;
  // if (z<=(normalDepth*alphaCoeff)) {
  if (z<=zo) {
    // shear zone
    return ((1.0-pow((1.0-z/zo),((nHB+1.0)/nHB)))*normalMaxVel);
  }
  else {
    // plug zone
    return (normalMaxVel);
  }
}

// double hbProfileBoundary (Point point, scalar s, double z, double ht)
// {
//   // general normal-flow profile of Herschel-Bulkley fluids
//   double zo;
//   zo = ht*(1.-bParam);
//   // zo = normalDepth*alphaCoeff;
//   // if (z<=(normalDepth*alphaCoeff)) {
//   if (z<=zo) {
//     // shear zone
//     return ((1.0-pow((1.0-z/zo),((nHB+1.0)/nHB)))*normalMaxVel);
//   }
//   else {
//     // plug zone
//     return (normalMaxVel);
//   }
// }

char s[80];
scalar depthGrad[], uAve[], yieldSurf[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = grav;
  init_grid(1 << (INITLEVEL));
  nl = 64;
  nu = 1.; // dummy

  CFL = 0.399;

//   periodic (right);

  run();
}

/**
## Initialization  */
event init (i = 0) {
//   double distDepth = normalDepth*(1.0+disAmp*sin(2. * pi * t / distPeriod));
  h[left] = t<=distPeriod/2.0 ? dirichlet(normalDepth*(1.0+disAmp*sin(2. * pi * t / distPeriod))) : dirichlet(normalDepth);
  h[right] = neumann(0.);

  for (vector u in ul) {
    for (int l = 0; l < nl; l++) {
      double z1 = 0.0;
      // FIXME: hardcore bc spec
      z1 += layer[l]*normalDepth/2.;
      u.n[left] = dirichlet(hbProfile(z, normalDepth));
      z1 += layer[l]*normalDepth/2.;
    }
    u.n[right] = neumann(0.);
  }

  /**
  We initialize *h*. */
  foreach(){
    double z = zb[];
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = normalDepth;


    uAve[] = 0.0;

    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);

      for (int l = 0; l < nl; l++) {
        z += layer[l]*h[]/2.;
        u = ul[l];
        u.x[] = hbProfile(z, normalDepth);
        // TODO: close follow Gaussian test case
//         ul.n[left] = t<=distPeriod/2.0 ? dirichlet(hbProfileBoundary(point, _s, z, distDepth)) : dirichlet(hbProfileBoundary(point, _s, z, normalDepth));
//         u.n[right] = neumann(0);
        z += layer[l]*h[]/2.;

        uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];
      }

    }
}

  /**
## Output
  We print the elevation  */
event acceleration (i++) {
  foreach(){
    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    uAve[] = 0.0;

      for (int l = 0; l < nl; l++) {
        u = ul[l];
        u.x[] += dt*grav;

        uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];
      }

  }
    boundary ((scalar *){u});
}

/*
 * AMR here
 *
 */
event adapt1 (i++) {
  adapt_wavelet({h, depthGrad, uAve}, (double[]){normalDepth/225.0, 0.010, normalVel/225.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   adapt_wavelet({h, depthGrad}, (double[]){normalDepth/250.0, 0.01250}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  refine(x<=refLen && t<=(2.0*distPeriod) && level<MAXLEVEL);
}

event hmax(i+=25)
{
     double maxDepth = normalDepth;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     double minDepth = normalDepth;
     double minDepthLocX = 0.0;
	 double max2LocX = 0.0;

     FILE *fp2 = fopen("maxMinDepth", "a+");

     // FR: maxDepth, minDepth, H, waveLength
     foreach ()
     {
               if ( h[]>h[1] && h[]>h[-1] && x>maxDepthLocX && h[]>(normalDepth*(1.0+disAmp))){
                    maxDepth = h[];
                    maxDepthLocX = x;
                    maxDepthVel = u.x[];
               }
     }
	 foreach ()
	 {
// 		if (h[]<h[1] && h[]<h[-1] && h[-1]<h[-2] && h[1]<h[2] && x>minDepthLocX && x<maxDepthLocX ) {
       if (h[]<h[1] && h[]<h[-1] && h[-1]<h[-2] && h[1]<h[2] && x>minDepthLocX && h[]<(normalDepth*(1.0+disAmp/2.0)) ) {
                    minDepthLocX = x;
                    minDepth = h[];
       }

// 		else if (h[]>h[1] && h[]>h[-1] && h[1]>h[2] && h[-1]>h[-2] && h[]>1.1 && x<maxDepthLocX && x>max2LocX) {
//                     max2LocX = x;
//                }
	 }
//      fprintf(fp2, "%g %13.9f %13.9f %13.9f %13.9f %13.9f \n", t, maxDepthLocX, maxDepth, minDepth, (maxDepth-minDepth), (maxDepthLocX-max2LocX));
    fprintf(fp2, "%g %13.9f %13.9f %13.9f %13.9f \n", t, maxDepthLocX, maxDepth, minDepth, (maxDepth-minDepth));
    fclose(fp2);
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval){
  sprintf (s, "slice-%g.txt", t);
  FILE * fp1 = fopen (s, "w"); 
  foreach(serial){
    double zCoord = zb[];
    double uVertGrad = 0.0;
    double uVertGradPrevLayer = 0.0;
    double detectYS = 0.0;

    uAve[] = 0.0;

    for (int l = 0; l < nl; l++) {
      zCoord += layer[l]*h[]*0.50;
      u = ul[l];

      if (l>0 && l<(nl-1)) {
        // middle layers
        vector um = ul[l-1] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l+1]*h[]*0.50+layer[l]*h[]);
      }
      else if (l==(nl-1))
      {
        // top layer
        vector um = ul[l-1] ;
        vector up = ul[l] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l]*h[]*1.50);
      }
      else
      {
        // bottom layer
        vector um = ul[l] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l+1]*h[]*0.50+layer[l]*h[]*1.50);
      }

      if (uVertGrad<yieldSurfThre && detectYS<1.0)
      {
        detectYS = 2.0;
        yieldSurf[] = zCoord-(layer[l-1]*h[]*0.50+layer[l]*h[]*0.50)*fabs((yieldSurfThre-uVertGrad)/(uVertGradPrevLayer-uVertGrad));
      }

      fprintf (fp1, "%g %g %g %g \n", x, zCoord, u.x[], uVertGrad);

      zCoord += layer[l]*h[]*0.50;

      uVertGradPrevLayer = uVertGrad;

      uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];
    }
  }
  fclose(fp1);

  sprintf (s, "depth-%g.txt", t);
  FILE * fp2 = fopen (s, "w"); 
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g \n", x, h[], uAve[], yieldSurf[]);
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
// 	   "# jet colormap\n"
// 	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
// 	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
// 	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
//        "# rainbow colormap\n"
       "set palette rgb 33,13,10\n"
      //  "load 'turbo.pal'\n"
//        "# green-red-violet colormap\n"
//        "set palette rgb 3,11,6\n"
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
      fprintf (fp, "%g %g %g \n", x, z, u.x[]);
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
