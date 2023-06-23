// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "saint-venant-power-law.h"
#include "utils.h"

// #define FR (1.756)
#define n_coeff 0.30
#define mun 150.0
#define rho 2130.0
#define grav 9.81
#define so 0.10
#define ALPHA (1/pow(FR,2.0))

#define aspectRatio 1.00
#define initDepth 1.00
// #define normalDepth 0.38093709
// #define normalVel 3.3915176

#define velThreshold 5.e-5

// double alpha_coeff = 40.0;
// double lx = 2.878;
// #define disMag 0.20
// double betaCoeff = 0.0;
#define simTime 3500.0
// double distPeriod = 1.0;
// #define wavelength (3.0/pow(FR,2.0))
// #define distPeriod (2.0/(FR*FR))
#define OUTPUTTIME1 0.2
#define OUTPUTTIME2 2.00
#define OUTPUTTIME3 50.00
#define OUTPUTPERIOD 20.0
#define OUTPUTPERIOD2 200.0
#define OUTPUTPLOT1 1.00
#define PLOTPERIOD 25.0
#define OUTPUTPLOT2 25.00
#define DOMAINLENGTH (initDepth*110.0)

#define XFROPERIOD 4.0

#define initialStageTime 3.0
#define restrictDt 1.0e-4

#define inletRefLen (aspectRatio*2.50)

// #define INITDEPTHFUNC(amp, len, xCoord) (1.0*(1.0+amp*sin(2.0*M_PI*(xCoord-len/2.0)/len)))
// #define INLETDEPTHBC(amp, time, period) (1.0*(1.0+amp*sin(2.0*M_PI*(time)/period)))

#define MAXLEVEL 14
#define MINLEVEL 2
#define INITLEVEL 14

#define WARMUP (0.00)

// h[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? INLETDEPTHBC(disMag, (t-WARMUP), distPeriod) : 1.0 );
// u.n[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? pow(INLETDEPTHBC(disMag, (t-WARMUP), distPeriod), 0.5) : 1.0 );

u.n[right] = neumann(0.);
h[right] = neumann(0.);

int main()
{
     L0 = DOMAINLENGTH;
//      N = 2048;
     init_grid (1 << INITLEVEL);
     // alphaCoeff = ALPHA;
     G = grav*1.0/pow((1+so*so),0.50);
     betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
     // nCoeff = n_coeff;
     
     CFL = 0.399; // CFL number should be sufficiently small

     // theta = 1.205;

     // note that we are using the hump periodic configuration
//      periodic (right);

     run();
}

scalar re[], fr[], frGrad[];
scalar te[], pe[], ke[];
scalar depthGrad[];

event init(i = 0)
{
     foreach ()
     {
               zb[] = 0.0;
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               // h[] = 1.0;
//                h[] = (x>=wavelength/2.0 && x<=wavelength) ? INITDEPTHFUNC(disMag, wavelength, x) : 1.0;
               // h[] = 1.0;
               h[] = x<=(initDepth*aspectRatio) ? initDepth : 0.0;
//                u.x[] = (x>=wavelength/2.0 && x<=wavelength) ? pow(INITDEPTHFUNC(disMag, wavelength, x), 0.5) : 1.0;
               // u.x[] = 1.0;
               u.x[] = 0.0;

               depthGrad[] = fabs(h[-1]-h[1])/Delta;

               fr[] = h[]>dry ? u.x[]/(pow(G*h[], 0.50)) : 0.0;
     }

}

event dtMax(t = 0; t<=initialStageTime; t+=restrictDt)
{

}

event calcDepthGrad(i++)
{
     foreach()
     {
          depthGrad[] = fabs(h[-1]-h[1])/Delta;
          frGrad[] = fabs(fr[-1]-fr[1])/Delta;
     }
}

static double powerLawFriction(double u, double h, double n)
{
     double rhs;
     rhs = G*so - mun/rho*pow(((1.+2.*n)/n*u/h), n)/h;
     return rhs;
}

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
          // if (h[] > dry && u.x[]>velThreshold ) {
          //      uMed = u.x[] + dt * powerLawFriction(u.x[], h[], n_coeff);
          //      uMed = (3.0/4.0)*u.x[] + (1.0/4.0)*uMed + (1.0/4.0)*dt*powerLawFriction(uMed, h[], n_coeff);
          //      u.x[] = (1.0/3.0)*u.x[]+(2.0/3.0)*uMed+(2.0/3.0)*dt*powerLawFriction(uMed, h[], n_coeff);
          // }
          // else {
          //      u.x[] = 0.0;
          // }

          //linearized backeard Euler
          u.x[] = h[]>dry && u.x[]>velThreshold  ? (u.x[] + dt*G*so)/(1.0 + dt*mun/rho*pow((1.+2.*n_coeff)/n_coeff,n_coeff)*pow(u.x[],(n_coeff-1.0))/pow(h[],(n_coeff+1.0))) : 0.0;
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));

          re[] = rho*pow(u.x[], (2.0-n_coeff))*pow(h[], n_coeff)/mun;
          fr[] = h[]>dry ? u.x[]/(pow(G*h[], 0.50)) : 0.0;
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

// record max depth
event hmaxUmax (i+=50)
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
                maxDepthVel = u.x[];
           }

           if (u.x[]>maxVel){
                maxVel = u.x[];
                maxVelLocX = x;
                maxVelDepth = h[];
           }
     }
     stats s1 = statsf (re);
     stats s2 = statsf (fr);
     fprintf(fp2, "%.10g %.10g %.10g %.10g %.10g \n", t, maxDepth, maxDepthLocX, maxDepthVel, s1.max);
     fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, maxVel, maxVelLocX, maxVelDepth);
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
     FILE *fp3 = fopen("frontPos", "a+");
     foreach ()
     {
          xf = h[] > dry ?  max(xf,x) :  xf ;
          aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
          aveVel += (h[]>=dry ? Delta*u.x[] : 0.0);
     }
//      stats s2 = statsf (frontXPos);
//      xf = s2.max;
     fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf));
     fclose(fp3);
}

event hFront2 (t=XFROPERIOD; i+=80)
{
     double aveDepth = 0.0;
     double aveVel = 0.0;
     double xf = 0.0;
     FILE *fp3 = fopen("frontPos", "a+");
     foreach ()
     {
          xf = h[] > dry ?  max(xf,x) :  xf ;
          aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
          aveVel += (h[]>=dry ? Delta*u.x[] : 0.0);
     }
//      stats s2 = statsf (frontXPos);
//      xf = s2.max;
     fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf));
     fclose(fp3);
}

event energyContent(i+=50)
{
     // Vallis textbook P65.
     FILE *fp2 = fopen("energyContent", "a+");
     foreach ()
     {
          te[] = 0.50*rho*h[]*pow(u.x[], 2.0) + 0.50*rho*G*h[];
          ke[] = 0.50*rho*h[]*pow(u.x[], 2.0);
          pe[] = 0.50*rho*G*h[];
     }
     stats s1 = statsf (te);
     stats s2 = statsf (ke);
     stats s3 = statsf (pe);
     fprintf(fp2, "%g %g %g %g\n", t, s1.sum, s2.sum, s3.sum);
     fclose(fp2);
}

event slopeCalc(i+=50)
{
     double slopeVal = 0.0;
     double endXCoord = 0.0;
     double beginDepth = 0.0;
     double endDepth = 0.0;
     int count = 0;
     FILE *fp1 = fopen("slope", "a+");
     foreach()
     {
          if(count<1)
          {
               beginDepth = h[];
               endDepth = h[];
          }
          else if (h[1]-h[-1]<0 && h[2]-h[]<0)
          {
               endDepth = h[];
               endXCoord = x;
               break;
          }
          count += 1;
     }
     slopeVal = (endDepth-beginDepth)/(endXCoord-0.0);
     fprintf(fp1, "%g %g \n", t, slopeVal);
     fclose(fp1);
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

event gnuplotOutput1(t = 0; t <= PLOTPERIOD; t += OUTPUTPLOT1)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
}

event gnuplotOutput2(t = PLOTPERIOD; t += OUTPUTPLOT2)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
}

event output1(t = 0.0; t <= OUTPUTPERIOD; t += OUTPUTTIME1)
{
     char name[30];
     sprintf(name, "out-%g.txt", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g %g %g \n", x, h[], u.x[], re[], fr[]);
//      fprintf(fp, "\n");
     fclose(fp);
}

event output2(t = OUTPUTPERIOD; t <= OUTPUTPERIOD2; t += OUTPUTTIME2)
{
     char name[30];
     sprintf(name, "out-%g.txt", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g %g %g \n", x, h[], u.x[], re[], fr[]);
//      fprintf(fp, "\n");
     fclose(fp);
}

event output3(t = OUTPUTPERIOD2; t <= simTime; t += OUTPUTTIME3)
{
     char name[30];
     sprintf(name, "out-%g.txt", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g %g %g \n", x, h[], u.x[], re[], fr[]);
//      fprintf(fp, "\n");
     fclose(fp);
}

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

//TODO: add AMR features
// event adapt (i++) {
// }
// AMR features
event adapt1 (i++) {
     adapt_wavelet({h, depthGrad, frGrad}, (double[]){initDepth/325.0, 1.8e-5, 4.e-4}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
     refine( x<=inletRefLen && t<=(pow(initDepth/G,0.5)*5.0) && level < MAXLEVEL);
}
