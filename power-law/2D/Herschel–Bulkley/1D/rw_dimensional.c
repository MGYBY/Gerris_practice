// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
// #include "saint-venant-power-law.h"
#include "conservation.h"
#include "utils.h"

// here we use the dimensionless form of governing eqns
scalar h[], q[];
scalar * scalars = {h,q};
vector * vectors = NULL;

#define FR (1.756)
#define nCoeff (1.00)
#define grav 9.81
#define so 0.06
#define gp (grav*pow((1.0-so*so),0.50))
#define tauC 200.0
#define rhoFluid 2130.0
#define muN 3.195

#define normalDepth 0.29326308
#define normalVel 2.9757479

// double alpha_coeff = 40.0;
// double lx = 2.878;
#define disMag 0.20
// double betaCoeff = 0.0;
#define simTime 60.0
// double distPeriod = 1.0;
#define wavelength (3.0/pow(FR,2.0))
#define distPeriod (2.0/(FR*FR))
#define OUTPUTTIME 2.00
#define DOMAINLENGTH (2.0*normalDepth/so)

#define uThreshold (normalVel/1000.0)
#define dryThreshold (normalDepth/1000.0)

#define inletRefLen (DOMAINLENGTH/20.0)

// #define INITDEPTHFUNC(amp, len, xCoord) (1.0*(1.0+amp*sin(2.0*M_PI*(xCoord-len/2.0)/len)))
#define INLETDEPTHBC(amp, time, period) (1.0*(1.0+amp*sin(2.0*M_PI*(time)/period)))

#define MAXLEVEL 12
#define MINLEVEL 4
// #define MAXMAXLEVEL 16
#define UNIFORMLEVEL 9

#define WARMUP (0.00)

void flux (const double * s, double * f, double e[2])
{
     double h = s[0], q = s[1], u = q/h;
     double alphaCoeff = (gp*h*(1 + 2*nCoeff)*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*gp*h*pow(1 + nCoeff,2)*so*rhoFluid))/((2 + 3*nCoeff)*pow(nCoeff*tauC + gp*h*(1 + nCoeff)*so*rhoFluid,2));
     f[0] = q;
     f[1] = alphaCoeff*q*u + gp*h*h/2.;

     // min/max eigenvalues
     // double c = sqrt(G*h);
     e[0] = (gp*(1 + 2*nCoeff)*q*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*gp*h*pow(1 + nCoeff,2)*so*rhoFluid) -
     sqrt(gp*(h*pow(nCoeff,4)*pow(2 + 3*nCoeff,2)*pow(tauC,4) +
         4*gp*pow(h,2)*pow(nCoeff,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2)*so*pow(tauC,3)*rhoFluid +
         gp*pow(nCoeff,2)*(pow(q,2) + 2*(1 + nCoeff)*
             (3*gp*pow(h,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(6 + nCoeff*(15 + 14*nCoeff))*pow(q,2)))*pow(so,2)*
          pow(tauC,2)*pow(rhoFluid,2) + 2*pow(gp,2)*h*nCoeff*pow(1 + nCoeff,2)*
          (2*gp*pow(h,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(1 + 2*nCoeff)*(3 + 4*nCoeff)*pow(q,2))*pow(so,3)*tauC*
          pow(rhoFluid,3) + pow(gp,3)*pow(h,2)*pow(1 + nCoeff,4)*
          (gp*pow(h,3)*pow(2 + 3*nCoeff,2) + 2*nCoeff*(1 + 2*nCoeff)*pow(q,2))*pow(so,4)*pow(rhoFluid,4))))/
   ((2 + 3*nCoeff)*pow(nCoeff*tauC + gp*h*(1 + nCoeff)*so*rhoFluid,2)); // min

     e[1] = (gp*(1 + 2*nCoeff)*q*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*gp*h*pow(1 + nCoeff,2)*so*rhoFluid) +
     sqrt(gp*(h*pow(nCoeff,4)*pow(2 + 3*nCoeff,2)*pow(tauC,4) +
         4*gp*pow(h,2)*pow(nCoeff,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2)*so*pow(tauC,3)*rhoFluid +
         gp*pow(nCoeff,2)*(pow(q,2) + 2*(1 + nCoeff)*
             (3*gp*pow(h,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(6 + nCoeff*(15 + 14*nCoeff))*pow(q,2)))*pow(so,2)*
          pow(tauC,2)*pow(rhoFluid,2) + 2*pow(gp,2)*h*nCoeff*pow(1 + nCoeff,2)*
          (2*gp*pow(h,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(1 + 2*nCoeff)*(3 + 4*nCoeff)*pow(q,2))*pow(so,3)*tauC*
          pow(rhoFluid,3) + pow(gp,3)*pow(h,2)*pow(1 + nCoeff,4)*
          (gp*pow(h,3)*pow(2 + 3*nCoeff,2) + 2*nCoeff*(1 + 2*nCoeff)*pow(q,2))*pow(so,4)*pow(rhoFluid,4))))/
   ((2 + 3*nCoeff)*pow(nCoeff*tauC + gp*h*(1 + nCoeff)*so*rhoFluid,2)); // max
}

// h[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? INLETDEPTHBC(disMag, (t-WARMUP), distPeriod) : 1.0 );
// u.n[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? pow(INLETDEPTHBC(disMag, (t-WARMUP), distPeriod), 0.5) : 1.0 );

// u.n[right] = neumann(0.);
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
     
     CFL = 0.399; // CFL number should be sufficiently small

     // theta = 1.205;

     // note that we are using the hump periodic configuration
     periodic (right);

     run();
}

scalar te[], re[];
// scalar inletRef[];
scalar depthGrad[], inletRef[];
scalar u[];

event init(i = 0)
{
     foreach ()
     {
          // zb[] = 0.0;
          // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
          // h[] = 1.0;
          // h[] = (x>=wavelength/2.0 && x<=wavelength) ? INITDEPTHFUNC(disMag, wavelength, x) : 1.0;
          // h[] = 1.0;
          h[] = normalDepth*(1.0+disMag*sin(2.0*M_PI*x/DOMAINLENGTH));
          // u.x[] = (x>=wavelength/2.0 && x<=wavelength) ? pow(INITDEPTHFUNC(disMag, wavelength, x), 0.5) : 1.0;
          // u.x[] = 1.0;
          q[] = normalDepth*normalVel;
          u[] = q[]/h[];

          // inletRef[] = x<=inletRefLen ? 1.0 : 0.0;
          depthGrad[] = fabs(h[-1]-h[1])/Delta;
     }

}

event calcDepthGrad(i++)
{
     foreach()
     {
          depthGrad[] = fabs(h[-1]-h[1])/Delta;
          // inletRef[] = (t<=5 && x<=inletRefLen) ? 1.0 : 0.0;
          u[] = (h[]>dryThreshold && u[]>uThreshold) ? q[]/h[] : 0.0;
     }
}

static double hbFriction(double q, double h)
{
     double rhs;
     rhs = gp*so*h-(1./rhoFluid)*(tauC+muN*pow(((q*pow((rhoFluid*gp*so),2.0)*(nCoeff+1.)*(2.*nCoeff+1.))/((rhoFluid*gp*so*h-tauC)*(nCoeff*(nCoeff+1.)*rhoFluid*h*gp*so+nCoeff*nCoeff*tauC))), nCoeff));
     return rhs;
}

// bottom friction
event friction(i++)
{
     double qMed=0.0;
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
          if (h[] > dryThreshold && u[]>uThreshold) {
               qMed = q[] + dt * hbFriction(q[], h[]);
               qMed = (3.0/4.0)*q[] + (1.0/4.0)*qMed + (1.0/4.0)*dt*hbFriction(qMed, h[]);
               q[] = (1.0/3.0)*q[]+(2.0/3.0)*qMed+(2.0/3.0)*dt*hbFriction(qMed, h[]);
          }
          else {
               q[] = 0.0;
          }

          //linearized backeard Euler
          // q[] = (h[]>dryThreshold && u[]>uThreshold) ? (0.0) : 0.0;
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));

          // re[] = pow(u.x[], (2.0-n_coeff))*pow(h[], n_coeff);
     }
     boundary((scalar *){q}); // note that the input should be a list (at least for 1d)
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
                    maxDepthVel = u[];
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
          // fprintf(fp, "%g %g %g %g \n", x, h[], u.x[], re[]);
          fprintf(fp, "%g %g %g %g\n", x, h[], q[], u[]);
     fprintf(fp, "\n");
     fclose(fp);
}

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

//TODO: add AMR features
// event adapt (i++) {
// }
// AMR features
// event adapt1 (i++) {
//   astats s = adapt_wavelet({h, depthGrad}, (double[]){1/225.0, 0.020}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
// //      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
// }
