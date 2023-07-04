/* A collection of Riemann solvers for the Saint-Venant system 
* power-law model
 */

#define epsilon 1e-30

void hlle (double hm, double hp, double um, double up, double Delta,
	   double * fh, double * fq, double * dtmax)
{
  // Roe average
  double uhat = (sqrt(hm)*um+sqrt(hp)*up)/(sqrt(hm)+sqrt(hp));
  double cm = sqrt (G*hm+sq(um)*(betaCoeff*betaCoeff-betaCoeff)), cp = sqrt (G*hp+sq(up)*(betaCoeff*betaCoeff-betaCoeff));
  double chat = sqrt(G*(hp+hm)/2.0+sq(uhat)*(betaCoeff*betaCoeff-betaCoeff));
  double SL = min (betaCoeff*um - cm, betaCoeff*uhat - chat); SL = min(SL, betaCoeff*up - cp);
  double SR = max (betaCoeff*up + cp, betaCoeff*uhat + chat); SR = max(SR, betaCoeff*um + cm);

  if (0. <= SL) {
    *fh = um*hm;
    *fq = hm*(betaCoeff*um*um + G*hm/2.);
  }
  else if (0. >= SR) {
    *fh = up*hp;
    *fq = hp*(betaCoeff*up*up + G*hp/2.);
  }
  else {
    double fhm = um*hm;
    double fum = hm*(betaCoeff*um*um + G*hm/2.);
    double fhp = up*hp;
    double fup = hp*(betaCoeff*up*up + G*hp/2.);
    *fh = (SR*fhm - SL*fhp + SL*SR*(hp - hm))/(SR - SL);
    *fq = (SR*fum - SL*fup + SL*SR*(hp*up - hm*um))/(SR - SL);
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > epsilon) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}
