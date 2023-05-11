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
  double alphaVal = (2.0*nCoeff+1.0)/(3.0*nCoeff+2.0)*(2.0*pow((nCoeff+1.0),2.0)*h+nCoeff*(4.0*nCoeff+3.0)*cStar)/(pow((nCoeff+1.0),2.0)*h+2.0*nCoeff*(nCoeff+1.0)*cStar+pow((nCoeff*cStar),2.0)/h);
  return alphaVal;
}

void kurganov (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double l1P = ((1 + 2*nCoeff)*(hp + (cStar + hp)*nCoeff)*(2*hp*pow(1 + nCoeff,2) + cStar*nCoeff*(3 + 4*nCoeff))*(hp*up) +
     sqrt(G*hp*pow(2 + 3*nCoeff,2)*pow(hp + (cStar + hp)*nCoeff,6) +
       nCoeff*(1 + 2*nCoeff)*pow(hp + (cStar + hp)*nCoeff,2)*
        (2*pow(hp,2)*pow(1 + nCoeff,4) + 2*cStar*hp*nCoeff*pow(1 + nCoeff,2)*(3 + 4*nCoeff) +
          pow(cStar,2)*nCoeff*(1 + 2*nCoeff*(5 + nCoeff*(11 + 7*nCoeff))))*pow((hp*up),2)))/
   ((2 + 3*nCoeff)*pow(hp + (cStar + hp)*nCoeff,3));

  double l1M = ((1 + 2*nCoeff)*(hm + (cStar + hm)*nCoeff)*(2*hm*pow(1 + nCoeff,2) + cStar*nCoeff*(3 + 4*nCoeff))*(hm*um) +
     sqrt(G*hm*pow(2 + 3*nCoeff,2)*pow(hm + (cStar + hm)*nCoeff,6) +
       nCoeff*(1 + 2*nCoeff)*pow(hm + (cStar + hm)*nCoeff,2)*
        (2*pow(hm,2)*pow(1 + nCoeff,4) + 2*cStar*hm*nCoeff*pow(1 + nCoeff,2)*(3 + 4*nCoeff) +
          pow(cStar,2)*nCoeff*(1 + 2*nCoeff*(5 + nCoeff*(11 + 7*nCoeff))))*pow((hm*um),2)))/
   ((2 + 3*nCoeff)*pow(hm + (cStar + hm)*nCoeff,3));

  double l2P = ((1 + 2*nCoeff)*(hp + (cStar + hp)*nCoeff)*(2*hp*pow(1 + nCoeff,2) + cStar*nCoeff*(3 + 4*nCoeff))*(hp*up) -
     sqrt(G*hp*pow(2 + 3*nCoeff,2)*pow(hp + (cStar + hp)*nCoeff,6) +
       nCoeff*(1 + 2*nCoeff)*pow(hp + (cStar + hp)*nCoeff,2)*
        (2*pow(hp,2)*pow(1 + nCoeff,4) + 2*cStar*hp*nCoeff*pow(1 + nCoeff,2)*(3 + 4*nCoeff) +
          pow(cStar,2)*nCoeff*(1 + 2*nCoeff*(5 + nCoeff*(11 + 7*nCoeff))))*pow((hp*up),2)))/
   ((2 + 3*nCoeff)*pow(hp + (cStar + hp)*nCoeff,3));

  double l2M = ((1 + 2*nCoeff)*(hm + (cStar + hm)*nCoeff)*(2*hm*pow(1 + nCoeff,2) + cStar*nCoeff*(3 + 4*nCoeff))*(hm*um) -
     sqrt(G*hm*pow(2 + 3*nCoeff,2)*pow(hm + (cStar + hm)*nCoeff,6) +
       nCoeff*(1 + 2*nCoeff)*pow(hm + (cStar + hm)*nCoeff,2)*
        (2*pow(hm,2)*pow(1 + nCoeff,4) + 2*cStar*hm*nCoeff*pow(1 + nCoeff,2)*(3 + 4*nCoeff) +
          pow(cStar,2)*nCoeff*(1 + 2*nCoeff*(5 + nCoeff*(11 + 7*nCoeff))))*pow((hm*um),2)))/
   ((2 + 3*nCoeff)*pow(hm + (cStar + hm)*nCoeff,3));

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
