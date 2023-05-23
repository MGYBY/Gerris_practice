#include <complex.h>    /* Standard Library of Complex Numbers */
// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
// #include "saint-venant-power-law.h"
// #include "conservation.h"
#include "./myConservation.h"
#include "utils.h"

// here we use the dimensionless form of governing eqns
scalar h[], q[], up[];
scalar * scalars = {h,q,up};
vector * vectors = NULL;

#define alphaCoeffVal 0.30
#define betaCoeffVal 27.0

#define normalDepth 1.00
#define normalUp (0.5*pow((1.0-alphaCoeff),2.0))
#define normalQ (normalUp*(normalDepth-1.0/3.0*(1.0-alphaCoeff)))

#define disMag 0.05
#define simTime 100.0
// #define wavelength (3.0/pow(FR,2.0))
// #define distPeriod (2.0/(FR*FR))
#define OUTPUTTIME 2.00
#define DOMAINLENGTH (2*3.14159265/1.20)

// #define uThreshold ((0.5*pow((1.-alphaCoeff),2.)*(2.0/3.0-alphaCoeff/3.0))/1000.0)
// #define dryThreshold (normalDepth/1000.0)

// #define inletRefLen (DOMAINLENGTH/25.0)

// #define INITDEPTHFUNC(amp, len, xCoord) (1.0*(1.0+amp*sin(2.0*M_PI*(xCoord-len/2.0)/len)))
// #define INLETDEPTHBC(amp, time, period, H) (H*(1.0+amp*sin(2.0*M_PI*(time)/period)))

#define MAXLEVEL 12
#define MINLEVEL 4
// #define MAXMAXLEVEL 16
#define UNIFORMLEVEL 9

#define WARMUP (0.00)

void flux (const double * s, double * f, double e[2], double betaCoeff)
{
     double h = s[0], q = s[1], up = s[2];
     f[0] = q;
     f[1] = (h-7.0/15.0*(3.0*(h - q/up)))*pow(up,2.0)+pow(h,2.0)/(2.0*betaCoeff);
     f[2] = pow(up,2.0)/2.0+h/betaCoeff;

     double complex eig1 = (0.2*(4.*up*betaCoeff + (1.2599210498948732*betaCoeff*(25.*h + 3.*cpow(up,2)*betaCoeff))/
        cpow(4725.*q*cpow(betaCoeff,2) - 3375.*h*up*cpow(betaCoeff,2) - 54.*cpow(up,3)*cpow(betaCoeff,3) +
          csqrt(729.*cpow(betaCoeff,4)*cpow(-175.*q + 125.*h*up + 2.*cpow(up,3)*betaCoeff,2) +
            4.*cpow(-75.*h*betaCoeff - 9.*cpow(up,2)*cpow(betaCoeff,2),3)),0.3333333333333333) +
       0.26456684199469993*cpow(4725.*q*cpow(betaCoeff,2) - 3375.*h*up*cpow(betaCoeff,2) - 54.*cpow(up,3)*cpow(betaCoeff,3) +
          csqrt(729.*cpow(betaCoeff,4)*cpow(-175.*q + 125.*h*up + 2.*cpow(up,3)*betaCoeff,2) +
            4.*cpow(-75.*h*betaCoeff - 9.*cpow(up,2)*cpow(betaCoeff,2),3)),0.3333333333333333)))/betaCoeff;

     double complex eig2 = 0.8*up - (CMPLX(0.12599210498948732, 0.21822472719434427)*(25.*h + 3.*cpow(up,2)*betaCoeff))/
    cpow(4725.*q*cpow(betaCoeff,2) - 3375.*h*up*cpow(betaCoeff,2) - 54.*cpow(up,3)*cpow(betaCoeff,3) +
      csqrt(729.*cpow(betaCoeff,4)*cpow(-175.*q + 125.*h*up + 2.*cpow(up,3)*betaCoeff,2) +
        4.*cpow(-75.*h*betaCoeff - 9.*cpow(up,2)*cpow(betaCoeff,2),3)),0.3333333333333333) -
   (CMPLX(0.026456684199469994, -0.04582432123328676)*
      cpow(4725.*q*cpow(betaCoeff,2) - 3375.*h*up*cpow(betaCoeff,2) - 54.*cpow(up,3)*cpow(betaCoeff,3) +
        csqrt(729.*cpow(betaCoeff,4)*cpow(-175.*q + 125.*h*up + 2.*cpow(up,3)*betaCoeff,2) +
          4.*cpow(-75.*h*betaCoeff - 9.*cpow(up,2)*cpow(betaCoeff,2),3)),0.3333333333333333))/betaCoeff;

     double complex eig3 = 0.8*up - (CMPLX(0.12599210498948732, -0.21822472719434427)*(25.*h + 3.*cpow(up,2)*betaCoeff))/
    cpow(4725.*q*cpow(betaCoeff,2) - 3375.*h*up*cpow(betaCoeff,2) - 54.*cpow(up,3)*cpow(betaCoeff,3) +
      csqrt(729.*cpow(betaCoeff,4)*cpow(-175.*q + 125.*h*up + 2.*cpow(up,3)*betaCoeff,2) +
        4.*cpow(-75.*h*betaCoeff - 9.*cpow(up,2)*cpow(betaCoeff,2),3)),0.3333333333333333) -
   (CMPLX(0.026456684199469994, 0.04582432123328676)*
      cpow(4725.*q*cpow(betaCoeff,2) - 3375.*h*up*cpow(betaCoeff,2) - 54.*cpow(up,3)*cpow(betaCoeff,3) +
        csqrt(729.*cpow(betaCoeff,4)*cpow(-175.*q + 125.*h*up + 2.*cpow(up,3)*betaCoeff,2) +
          4.*cpow(-75.*h*betaCoeff - 9.*cpow(up,2)*cpow(betaCoeff,2),3)),0.3333333333333333))/betaCoeff;

     // min/max eigenvalues
     // double c = sqrt(G*h);
     e[0] = min(creal(eig1), min(creal(eig2), creal(eig3))); // min

     e[1] = max(creal(eig1), max(creal(eig2), creal(eig3))); // max
}

// h[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? INLETDEPTHBC(disMag, (t-WARMUP), distPeriod, normalDepth) : normalDepth );
// q[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? FR*pow(gp*pow(INLETDEPTHBC(disMag, (t-WARMUP), distPeriod, normalDepth), 3.0), 0.5) : (normalDepth*normalVel) );

// q[right] = neumann(0.);
// h[right] = neumann(0.);

int main()
{
     L0 = DOMAINLENGTH;
//      N = 2048;
     init_grid (1 << UNIFORMLEVEL);
     // alphaCoeff = ALPHA;
     // G = gp;
     // betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
     // nCoeff = n_coeff;

     betaCoeff = betaCoeffVal;
     alphaCoeff = alphaCoeffVal;
     
     CFL = 0.399; // CFL number should be sufficiently small

     theta = 1.20;

     // note that we are using the hump periodic configuration
     periodic (right);

     run();
}

scalar te[], re[];
// scalar inletRef[];
scalar depthGrad[], inletRef[];
// scalar u[];

event init(i = 0)
{
     foreach ()
     {
          // zb[] = 0.0;
          // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
          // h[] = 1.0;
          // h[] = (x>=wavelength/2.0 && x<=wavelength) ? INITDEPTHFUNC(disMag, wavelength, x) : 1.0;
          h[] = normalDepth*(1.+disMag*sin(2.0*M_PI*x/DOMAINLENGTH));
          // u.x[] = (x>=wavelength/2.0 && x<=wavelength) ? pow(INITDEPTHFUNC(disMag, wavelength, x), 0.5) : 1.0;
          // u.x[] = 1.0;
          q[] = normalQ;
          // u[] = q[]/h[];
          up[] = normalUp;

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
          // u[] = (h[]>dryThreshold && u[]>uThreshold) ? q[]/h[] : 0.0;
          // re[] = 8.0*rhoFluid*pow(u[],2.0)/(tauC+muN*pow((2.0*u[]/h[]),nCoeff));
     }

     fprintf(fp1, "%d \n", i);

     fclose(fp1);
}

static double signNum(double x) {
     if (x > 0.0) {
          return 1.;
     }
     else if (x < 0.0) {
          return -1.;
     }
     else {
          return 0.0;
     }
}

static double hbFriction1(double h, double q, double up)
{
     double rhs;
     rhs = 1./betaCoeff*(h-alphaCoeff*signNum(up)-2.*up/(3.*(h - q/up)));
     return rhs;
}

static double hbFriction2(double h, double q, double up)
{
     double rhs;
     rhs = 1./betaCoeff*(1.-alphaCoeff*signNum(up)/(h-3.*(h - q/up)));
     return rhs;
}

// bottom friction
event friction(i++)
{
     double qMed=0.0;
     double upMed=0.0;
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
          // wet bed && not too slow && ho<h
          if (h[] > (normalDepth/1000.0) && up[]>(normalUp/1000.0) && (3.0*(h[]-q[]/up[]))<h[]) {
               // Fixme: modify to exact formulation of RK3TVD
               qMed = q[] + dt * hbFriction1(h[], q[], up[]);
               qMed = (3.0/4.0)*q[] + (1.0/4.0)*qMed + (1.0/4.0)*dt*hbFriction1(h[], qMed, up[]);
               q[] = (1.0/3.0)*q[]+(2.0/3.0)*qMed+(2.0/3.0)*dt*hbFriction1(h[], qMed, up[]);

               upMed = up[] + dt * hbFriction2(h[], q[], up[]);
               upMed = (3.0/4.0)*up[] + (1.0/4.0)*upMed + (1.0/4.0)*dt*hbFriction2(h[], q[], upMed);
               up[] = (1.0/3.0)*up[]+(2.0/3.0)*upMed+(2.0/3.0)*dt*hbFriction2(h[], q[], upMed);
          }
          else {
               q[] = 0.0;
               up[] = 0.0;
               qMed = 0.0;
               upMed = 0.0;
          }

          //linearized backeard Euler
          // q[] = (h[]>dryThreshold && u[]>uThreshold) ? ((q[]+dt*(gp*h[]*so+(1./rhoFluid)*(-1.*tauC)))/(1.+dt*muN/rhoFluid*pow(q[],(nCoeff-1.))*pow((pow((rhoFluid*gp*so),2.0)*(nCoeff+1.)*(2.*nCoeff+1.))/((rhoFluid*gp*so*h[]-tauC)*(nCoeff*(nCoeff+1.)*rhoFluid*h[]*gp*so+nCoeff*nCoeff*tauC)),nCoeff))) : 0.0;
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));
          
     }
     boundary((scalar *){q, up}); // note that the input should be a list (at least for 1d)
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
               // maxDepthVel = u[];
          }
     }
     // stats s1 = statsf (re);
     // fprintf(fp2, "%g %g %g %g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel, s1.max);
     fprintf(fp2, "%g %g %g \n", t, maxDepth, maxDepthLocX);
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
          // fprintf(fp, "%g %g %g %g %g \n", x, h[], q[], u[], re[]);
          fprintf(fp, "%g %g %g %g\n", x, h[], q[], up[]);
     fprintf(fp, "\n");
     fclose(fp);
}

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

//TODO: add AMR features
// event adapt1 (i++) {
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){(normalDepth/250.0), 0.050}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
// //      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

//      refine(x<=inletRefLen && level < MAXLEVEL);
// }
