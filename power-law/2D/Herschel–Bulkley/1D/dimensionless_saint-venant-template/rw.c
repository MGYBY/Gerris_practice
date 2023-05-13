// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "saint-venant-hb.h"
#include "utils.h"

#define FR (0.50)
#define cStarVal 0.40
#define nCoeffVal 1.0

// double alpha_coeff = 40.0;
// double lx = 2.878;
#define disMag 0.20
// double betaCoeff = 0.0;
#define simTime 64.0
// double distPeriod = 1.0;
#define wavelength (3.0/pow(FR,2.0))
#define distPeriod (2.0/(FR*FR))
#define OUTPUTTIME 2.00
#define DOMAINLENGTH (200/(FR*FR))

#define uThreshold (1.0/1000.0)
#define dryThreshold (1.0e-5)

#define inletRefLen (distPeriod*5.0)

// #define INITDEPTHFUNC(amp, len, xCoord) (1.0*(1.0+amp*sin(2.0*M_PI*(xCoord-len/2.0)/len)))
#define INLETDEPTHBC(amp, time, period, H) (H*(1.0+amp*sin(2.0*M_PI*(time)/period)))

#define MAXLEVEL 16
#define MINLEVEL 4
// #define MAXMAXLEVEL 16
#define UNIFORMLEVEL 10

#define WARMUP (0.00)


h[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? INLETDEPTHBC(disMag, (t-WARMUP), distPeriod, 1.0) : (1.0) );
u.n[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? pow(INLETDEPTHBC(disMag, (t-WARMUP), distPeriod, 1.0), 0.5) : (1.0) );

u.n[right] = neumann(0.);
h[right] = neumann(0.);

int main()
{
     L0 = DOMAINLENGTH;
//      N = 2048;
     init_grid (1 << UNIFORMLEVEL);
     // alphaCoeff = ALPHA;
     G = 1.0/(FR*FR);
     cStar = cStarVal;
     nCoeff = nCoeffVal;
     
     CFL = 0.39; // CFL number should be sufficiently small

     // theta = 1.205;

     // note that we are using the hump periodic configuration
     // periodic (right);

     run();
}

scalar te[], re[];
scalar depthGrad[], inletRef[];

event init(i = 0)
{
     foreach ()
     {
          zb[] = 0.0;
          // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
          // h[] = 1.0;
          // h[] = (x>=wavelength/2.0 && x<=wavelength) ? INITDEPTHFUNC(disMag, wavelength, x) : 1.0;
          h[] = 1.0;
          // h[] = normalDepth*(1.0+disMag*sin(2.0*M_PI*x/DOMAINLENGTH));
          // u.x[] = (x>=wavelength/2.0 && x<=wavelength) ? pow(INITDEPTHFUNC(disMag, wavelength, x), 0.5) : 1.0;
          u.x[] = 1.0;
          // q[] = normalDepth*normalVel;

          // inletRef[] = x<=inletRefLen ? 1.0 : 0.0;
          depthGrad[] = fabs(h[-1]-h[1])/Delta;
     }

}

event calcDepthGrad(i++)
{
     FILE *fp1 = fopen("debugTimeStp", "a+");

     foreach()
     {
          depthGrad[] = fabs(h[-1]-h[1])/Delta;
          // inletRef[] = (t<=5 && x<=inletR) ? 1.0 : 0.0;
          // re[] = 8.0*rhoFluid*pow(u.x[],2.0)/(tauC+muN*pow((2.0*u.x[]/h[]),nCoeff));
     }

     fprintf(fp1, "%d \n", i);

     fclose(fp1);
}

// static double hbFriction(double u, double h)
// {
//      double rhs;
//      rhs = G*so-(1./rhoFluid)*(tauC/h+muN*pow(h,(nCoeff-1.))*pow(((u*pow((rhoFluid*G*so),2.0)*(nCoeff+1.)*(2.*nCoeff+1.))/((rhoFluid*G*so*h-tauC)*(nCoeff*(nCoeff+1.)*rhoFluid*h*G*so+nCoeff*nCoeff*tauC))), nCoeff));
//      return rhs;
// }

// bottom friction
event friction(i++)
{
     // double uMed=0.0;
     foreach ()
     {
          // rk2tvd
          // if (h[] > dry) {
          //      uMed = u.x[] + dt * powerLawFriction(u.x[], h[], n_coeff);
          //      u.x[] = 0.5*u.x[]+0.5*uMed+0.5*dt*powerLawFriction(uMed, h[], n_coeff);
          // }
          // else {
          //      u.x[] = 0.0;
          // }

          // rk3tvd
          // if (h[] > dryThreshold && u.x[]>uThreshold) {
          //      uMed = u.x[] + dt * hbFriction(u.x[], h[]);
          //      if (uMed>uThreshold)
          //      {
          //           uMed = (3.0/4.0)*u.x[] + (1.0/4.0)*uMed + (1.0/4.0)*dt*hbFriction(uMed, h[]);
          //           if (uMed>uThreshold)
          //           {
          //                u.x[] = (1.0/3.0)*u.x[]+(2.0/3.0)*uMed+(2.0/3.0)*dt*hbFriction(uMed, h[]);
          //           }
          //           else
          //           {
          //                u.x[] = 0.0;
          //           }
          //      }
          //      else {
          //           u.x[] = 0.0;
          //      }
          // }
          // else {
          //      u.x[] = 0.0;
          // }

          //linearized backward Euler
          u.x[] = ((u.x[]>uThreshold)&&(h[]>dryThreshold)&&(h[]>cStar)) ? ((u.x[]+dt*(1.0-cStar/h[]))/(1.0+dt*(1.-cStar)*pow((u.x[]*h[]),(nCoeff-1.))*pow((((1.-cStar)*(nCoeff+1.+nCoeff*cStar))/((h[]-cStar)*((nCoeff+1.)*h[]+nCoeff*nCoeff*cStar))),nCoeff))) : 0.0;
          
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));
          
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

// record max depth
event hmax(i+=20)
{
     double maxDepth = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     FILE *fp2 = fopen("maxDepth", "a+");
     foreach ()
     {
          if (h[]>maxDepth){
               maxDepth = h[];
               maxDepthLocX = x;
               maxDepthVel = u.x[];
          }
     }
     // stats s1 = statsf (re);
     // fprintf(fp2, "%g %g %g %g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel, s1.max);
     fprintf(fp2, "%g %g %g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel);
     fclose(fp2);
}

// event energyContent(i+=10)
// {
//      FILE *fp2 = fopen("totalEnergy", "a+");
//      foreach ()
//      {
//           te[] = 0.50*u.x[]*u.x[]+1.0/(pow(FR, 2.0))*h[];
//      }
//      stats s2 = statsf (te);
//      fprintf(fp2, "%g %g \n", t, s2.sum);
//      fclose(fp2);
// }

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

event gnuplot(t = 0; t <= simTime; t += OUTPUTTIME)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
     // fprintf(fp,
     //         "set term pngcairo enhanced size 800,600 font \",10\"\n"
     //         "set output 't%.0f.png'\n"
     //         "set title 't = %.2f'\n"
     //         "set xrange [0:40]\n"
     //         "plot u 1:2 w l t\n",
     //         t, t);
     // fprintf(fp, "\n");
     // foreach ()
     //      fprintf(fp, "%g %g\n", x, h[]);
     // fprintf(fp, "e\n\n");
     // fflush(fp);
     // fprintf(stderr, "%.3f %.3f\n", t, statsf(h).max); // uncomment if needed
}

event output(t = 0; t <= simTime; t += OUTPUTTIME)
{
     char name[80];
     sprintf(name, "out-%g.txt", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g \n", x, h[], u.x[]);
          // fprintf(fp, "%g %g %g %g\n", x, h[], q[], u[]);
     fprintf(fp, "\n");
     fclose(fp);
}

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

//TODO: add AMR features
event adapt1 (i++) {
     astats s = adapt_wavelet({h, depthGrad}, (double[]){(1.0/275.0), 0.0130}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
     fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

     refine(x<=inletRefLen && t<=(distPeriod*2.70) && level < MAXLEVEL);
}
