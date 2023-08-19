/* A collection of Riemann solvers for the Saint-Venant system 
* power-law model
 */

#define epsilon 1e-30

static double mmo (double a, double b)
{
  if (a>0 && b>0) {
    if (a>b) return b;
    else return a;
  }
  else if (a<0 && b<0) {
    if (a>b) return a;
    else return b;
  }
  else {
    return 0.0;
  }
}

void kurganov (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double cp = sqrt(G*hp+sq(up)*(betaCoeff*betaCoeff-betaCoeff)), cm = sqrt(G*hm+sq(um)*(betaCoeff*betaCoeff-betaCoeff));
  double ap = max(betaCoeff*up + cp, betaCoeff*um + cm); ap = max(ap, 0.);
  double am = min(betaCoeff*up - cp, betaCoeff*um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);
  if (a > epsilon) {
    *fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fq = (ap*(betaCoeff*qm*um + G*sq(hm)/2.) - am*(betaCoeff*qp*up + G*sq(hp)/2.) +
	    ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = 0.;
}

void kurganovRH (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double tinyParam = 1.0e-10;
  double cp = sqrt(G*hp+sq(up)*(betaCoeff*betaCoeff-betaCoeff)), cm = sqrt(G*hm+sq(um)*(betaCoeff*betaCoeff-betaCoeff));
  double ap = max(betaCoeff*up + cp, betaCoeff*um + cm); ap = max(ap, 0.);
  double am = min(betaCoeff*up - cp, betaCoeff*um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double deltaUEps1 = (hp-hm)>0 ? max((hp-hm), tinyParam) : min((hp-hm), -1.0*tinyParam);
  double deltaUEps2 = (qp-qm)>0 ? max((qp-qm), tinyParam) : min((qp-qm), -1.0*tinyParam);
  double shat1 = (2.0*(qp-qm))/((hp-hm)+deltaUEps1);
  double shat2 = (2.0*((betaCoeff*qp*up + G*pow(hp ,2.0)/2.)-(betaCoeff*qm*um + G*pow(hm ,2.0)/2.)))/((qp-qm)+deltaUEps2);
  double smax = max(shat1, shat2), smin = min(shat1, shat2);
  if (smax>tinyParam)
  {
    ap = min(ap, smax);
    am = max(am, -1.0*smax);
  }
  else if (smin<(-1.0*tinyParam))
  {
    ap = min(ap, -1.0*smin);
    am = max(am, smin);
  }
  double a = max(ap, -am);
  if (a > epsilon) {
    double wint = (ap*hp-am*hm-(qp-qm))/(ap-am);
    double qCorr = mmo((hp-wint)/(ap-am), (wint-hm)/(ap-am));
    *fh = (ap*qm - am*qp)/(ap - am) + (ap*am*((hp - hm)/(ap - am)-qCorr)); // (4.5) of [1]
    wint = (ap*qp-am*qm-((betaCoeff*qp*up + G*pow(hp ,2.0)/2.)-(betaCoeff*qm*um + G*pow(hm ,2.0)/2.)))/(ap-am);
    qCorr = mmo((qp-wint)/(ap-am), (wint-qm)/(ap-am));
    *fq = (ap*(betaCoeff*qm*um + G*pow(hp ,2.0)/2.) - am*(betaCoeff*qp*up + G*pow(hp ,2.0)/2.))/(ap - am) + (ap*am*((qp - qm)/(ap - am)-qCorr));

    //*fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    //*fq = (ap*(1.0*qm*um + G*sq(hm)/2.) - am*(1.0*qp*up + G*sq(hp)/2.) +
	    //ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = 0.;
}
