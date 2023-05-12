/* A collection of Riemann solvers for the Saint-Venant system 
 *
 * References:
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */

#define epsilon 1e-30

static double alphaCoeff (double h)
{
  // double alphaVal = (G*h*(1. + 2.*nCoeff)*so*rhoFluid*(nCoeff*(3. + 4.*nCoeff)*tauC + 2.*G*h*pow(1 + nCoeff,2)*so*rhoFluid))/((2. + 3.*nCoeff)*pow(nCoeff*tauC + G*h*(1. + nCoeff)*so*rhoFluid,2.));
  double alphaVal = ((2.*nCoeff+1.)/(3.*nCoeff+2.))*((2.*pow((nCoeff+1.),2.)*rhoFluid*G*so*h+tauC*nCoeff*(4.*nCoeff+3.))/(pow((nCoeff+1.),2.)*rhoFluid*G*so*h+2.*nCoeff*(nCoeff+1.)+2.*nCoeff*(nCoeff+1.)*tauC+(pow((nCoeff*tauC),2.))/(rhoFluid*G*so*h)));
  return alphaVal;
}

void kurganov (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double l1P = (G*(1 + 2*nCoeff)*(hp*up)*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*G*hp*pow(1 + nCoeff,2)*so*rhoFluid) +
     sqrt(G*(hp*pow(nCoeff,4)*pow(2 + 3*nCoeff,2)*pow(tauC,4) +
         4*G*pow(hp,2)*pow(nCoeff,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2)*so*pow(tauC,3)*rhoFluid +
         G*pow(nCoeff,2)*(pow((hp*up),2) + 2*(1 + nCoeff)*
             (3*G*pow(hp,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(6 + nCoeff*(15 + 14*nCoeff))*pow((hp*up),2)))*pow(so,2)*
          pow(tauC,2)*pow(rhoFluid,2) + 2*pow(G,2)*hp*nCoeff*pow(1 + nCoeff,2)*
          (2*G*pow(hp,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(1 + 2*nCoeff)*(3 + 4*nCoeff)*pow((hp*up),2))*pow(so,3)*tauC*
          pow(rhoFluid,3) + pow(G,3)*pow(hp,2)*pow(1 + nCoeff,4)*
          (G*pow(hp,3)*pow(2 + 3*nCoeff,2) + 2*nCoeff*(1 + 2*nCoeff)*pow((hp*up),2))*pow(so,4)*pow(rhoFluid,4))))/
   ((2 + 3*nCoeff)*pow(nCoeff*tauC + G*hp*(1 + nCoeff)*so*rhoFluid,2));

  double l1M = (G*(1 + 2*nCoeff)*(hm*um)*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*G*hm*pow(1 + nCoeff,2)*so*rhoFluid) +
     sqrt(G*(hm*pow(nCoeff,4)*pow(2 + 3*nCoeff,2)*pow(tauC,4) +
         4*G*pow(hm,2)*pow(nCoeff,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2)*so*pow(tauC,3)*rhoFluid +
         G*pow(nCoeff,2)*(pow((hm*um),2) + 2*(1 + nCoeff)*
             (3*G*pow(hm,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(6 + nCoeff*(15 + 14*nCoeff))*pow((hm*um),2)))*pow(so,2)*
          pow(tauC,2)*pow(rhoFluid,2) + 2*pow(G,2)*hm*nCoeff*pow(1 + nCoeff,2)*
          (2*G*pow(hm,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(1 + 2*nCoeff)*(3 + 4*nCoeff)*pow((hm*um),2))*pow(so,3)*tauC*
          pow(rhoFluid,3) + pow(G,3)*pow(hm,2)*pow(1 + nCoeff,4)*
          (G*pow(hm,3)*pow(2 + 3*nCoeff,2) + 2*nCoeff*(1 + 2*nCoeff)*pow((hm*um),2))*pow(so,4)*pow(rhoFluid,4))))/
   ((2 + 3*nCoeff)*pow(nCoeff*tauC + G*hm*(1 + nCoeff)*so*rhoFluid,2));

  double l2P = (G*(1 + 2*nCoeff)*(hp*up)*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*G*hp*pow(1 + nCoeff,2)*so*rhoFluid) -
     sqrt(G*(hp*pow(nCoeff,4)*pow(2 + 3*nCoeff,2)*pow(tauC,4) +
         4*G*pow(hp,2)*pow(nCoeff,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2)*so*pow(tauC,3)*rhoFluid +
         G*pow(nCoeff,2)*(pow((hp*up),2) + 2*(1 + nCoeff)*
             (3*G*pow(hp,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(6 + nCoeff*(15 + 14*nCoeff))*pow((hp*up),2)))*pow(so,2)*
          pow(tauC,2)*pow(rhoFluid,2) + 2*pow(G,2)*hp*nCoeff*pow(1 + nCoeff,2)*
          (2*G*pow(hp,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(1 + 2*nCoeff)*(3 + 4*nCoeff)*pow((hp*up),2))*pow(so,3)*tauC*
          pow(rhoFluid,3) + pow(G,3)*pow(hp,2)*pow(1 + nCoeff,4)*
          (G*pow(hp,3)*pow(2 + 3*nCoeff,2) + 2*nCoeff*(1 + 2*nCoeff)*pow((hp*up),2))*pow(so,4)*pow(rhoFluid,4))))/
   ((2 + 3*nCoeff)*pow(nCoeff*tauC + G*hp*(1 + nCoeff)*so*rhoFluid,2));

  double l2M = (G*(1 + 2*nCoeff)*(hm*um)*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*G*hm*pow(1 + nCoeff,2)*so*rhoFluid) -
     sqrt(G*(hm*pow(nCoeff,4)*pow(2 + 3*nCoeff,2)*pow(tauC,4) +
         4*G*pow(hm,2)*pow(nCoeff,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2)*so*pow(tauC,3)*rhoFluid +
         G*pow(nCoeff,2)*(pow((hm*um),2) + 2*(1 + nCoeff)*
             (3*G*pow(hm,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(6 + nCoeff*(15 + 14*nCoeff))*pow((hm*um),2)))*pow(so,2)*
          pow(tauC,2)*pow(rhoFluid,2) + 2*pow(G,2)*hm*nCoeff*pow(1 + nCoeff,2)*
          (2*G*pow(hm,3)*(1 + nCoeff)*pow(2 + 3*nCoeff,2) + nCoeff*(1 + 2*nCoeff)*(3 + 4*nCoeff)*pow((hm*um),2))*pow(so,3)*tauC*
          pow(rhoFluid,3) + pow(G,3)*pow(hm,2)*pow(1 + nCoeff,4)*
          (G*pow(hm,3)*pow(2 + 3*nCoeff,2) + 2*nCoeff*(1 + 2*nCoeff)*pow((hm*um),2))*pow(so,4)*pow(rhoFluid,4))))/
   ((2 + 3*nCoeff)*pow(nCoeff*tauC + G*hm*(1 + nCoeff)*so*rhoFluid,2));

  // double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  // double am = min(up - cp, um - cm); am = min(am, 0.);
  double ap = max(l1P, l1M); ap = max(ap, 0.);
  double am = min(l2P, l2M); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);

  // double alphaCoeff = (G*h*(1 + 2*nCoeff)*so*rhoFluid*(nCoeff*(3 + 4*nCoeff)*tauC + 2*G*h*pow(1 + nCoeff,2)*so*rhoFluid))/((2 + 3*nCoeff)*pow(nCoeff*tauC + G*h*(1 + nCoeff)*so*rhoFluid,2));

  if (a > epsilon) {
    *fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fq = (ap*(alphaCoeff(hm)*qm*um + G*sq(hm)/2.) - am*(alphaCoeff(hp)*qp*up + G*sq(hp)/2.) + 
	    ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = 0.;
}
