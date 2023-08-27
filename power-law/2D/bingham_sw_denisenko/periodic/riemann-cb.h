/* A collection of Riemann solvers for the Saint-Venant system 
 *
 * References:
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */

#define SQRT3 1.73205080756888
#define epsilon 1e-30

static double mmo (double a, double b)
{
  if (a>0 && b>0) {
    if (a>b) return b;
    if (b>a) return a;
  }
  else if (a<0 && b<0) {
    if (a>b) return a;
    if (b>a) return b;
  }
  else {
    return 0.0;
  }
}

void kurganov (double hm, double hp, double um, double up, double hiem, double hiep, double Delta,
	       double * fh, double * fq, double * fhie, double * dtmax)
{
  double phim = (2.*hiem-hm*(G*hm+pow(um,2.0)))/pow(hm,3.0), phip = (2.*hiep-hp*(G*hp+pow(up,2.0)))/pow(hp,3.0);
  double cp = sqrt(G*hp+3.*hp*hp*phip), cm = sqrt(G*hm+3.*hm*hm*phim);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);

  if (a > epsilon) {
    *fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fq = (ap*(qm*um + G*sq(hm)/2.+pow(hm,3.0)*phim) - am*(qp*up + G*sq(hp)/2.+pow(hp,3.0)*phip) +
	    ap*am*(qp - qm))/(ap - am);
    *fhie = (ap*(qm*(hiem/hm+(G*sq(hm)/2.+pow(hm,3.0)*phim)/hm)) - am*(qp*(hiep/hp+(G*sq(hp)/2.+pow(hp,3.0)*phip)/hp)) +
	    ap*am*(hiep - hiem))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = *fhie = 0.;
}

void kurganovSharp (double hm, double hp, double um, double up, double hiem, double hiep, double Delta,
	       double * fh, double * fq, double * fhie, double * dtmax)
{
  double phim = (2.*hiem-hm*(G*hm+pow(um,2.0)))/pow(hm,3.0), phip = (2.*hiep-hp*(G*hp+pow(up,2.0)))/pow(hp,3.0);
  double cp = sqrt(G*hp+3.*hp*hp*phip), cm = sqrt(G*hm+3.*hm*hm*phim);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);
  if (a > epsilon) {
    double wint = (ap*hp-am*hm-(qp-qm))/(ap-am);
    double qCorr = mmo((hp-wint)/(ap-am), (wint-hm)/(ap-am));
    *fh = (ap*qm - am*qp)/(ap - am) + (ap*am*((hp - hm)/(ap - am)-qCorr)); // (4.5) of [1]
    wint = (ap*qp-am*qm-((qp*up + G*pow(hp ,2.0)/2.+pow(hp,3.0)*phip)-(qm*um + G*pow(hm ,2.0)/2.+pow(hm,3.0)*phim)))/(ap-am);
    qCorr = mmo((qp-wint)/(ap-am), (wint-qm)/(ap-am));
    *fq = (ap*(qm*um + G*sq(hm)/2.+pow(hm,3.0)*phim) - am*(qp*up + G*sq(hp)/2.+pow(hp,3.0)*phip) )/(ap - am) + (ap*am*((qp - qm)/(ap - am)-qCorr));
    wint = (ap*hiep-am*hiem-((qp*(hiep/hp+(G*sq(hp)/2.+pow(hp,3.0)*phip)/hp))-(qm*(hiem/hm+(G*sq(hm)/2.+pow(hm,3.0)*phim)/hm))))/(ap-am);
    qCorr = mmo((hiep-wint)/(ap-am), (wint-hiem)/(ap-am));
    *fhie = (ap*(qm*(hiem/hm+(G*sq(hm)/2.+pow(hm,3.0)*phim)/hm)) - am*(qp*(hiep/hp+(G*sq(hp)/2.+pow(hp,3.0)*phip)/hp)))/(ap - am) + (ap*am*((hiep - hiem)/(ap - am)-qCorr));

    //*fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    //*fq = (ap*(1.0*qm*um + G*sq(hm)/2.) - am*(1.0*qp*up + G*sq(hp)/2.) +
	    //ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = *fhie = 0.;
}

// TODO: add the intuitive HLLC solver

// void hlle (double hm, double hp, double um, double up, double Delta,
// 	   double * fh, double * fq, double * dtmax)
// {
//   // Roe average
//   double uhat = (sqrt(hm)*um+sqrt(hp)*up)/(sqrt(hm)+sqrt(hp));
//   double cm = sqrt (G*hm+sq(um)*(1.0*1.0-1.0)), cp = sqrt (G*hp+sq(up)*(1.0*1.0-1.0));
//   double chat = sqrt(G*(hp+hm)/2.0+sq(uhat)*(1.0*1.0-1.0));
//   double SL = min (1.0*um - cm, 1.0*uhat - chat); SL = min(SL, 1.0*up - cp);
//   double SR = max (1.0*up + cp, 1.0*uhat + chat); SR = max(SR, 1.0*um + cm);
//
//   if (0. <= SL) {
//     *fh = um*hm;
//     *fq = hm*(1.0*um*um + G*hm/2.);
//   }
//   else if (0. >= SR) {
//     *fh = up*hp;
//     *fq = hp*(1.0*up*up + G*hp/2.);
//   }
//   else {
//     double fhm = um*hm;
//     double fum = hm*(1.0*um*um + G*hm/2.);
//     double fhp = up*hp;
//     double fup = hp*(1.0*up*up + G*hp/2.);
//     *fh = (SR*fhm - SL*fhp + SL*SR*(hp - hm))/(SR - SL);
//     *fq = (SR*fum - SL*fup + SL*SR*(hp*up - hm*um))/(SR - SL);
//   }
//
//   double a = max(fabs(SL), fabs(SR));
//   if (a > epsilon) {
//     double dt = CFL*Delta/a;
//     if (dt < *dtmax)
//       *dtmax = dt;
//   }
// }

// void kurganovRH (double hm, double hp, double um, double up, double Delta,
// 	       double * fh, double * fq, double * dtmax)
// {
//   double tinyParam = 1.0e-10;
//   double cp = sqrt(G*hp+sq(up)*(1.0*1.0-1.0)), cm = sqrt(G*hm+sq(um)*(1.0*1.0-1.0));
//   double ap = max(1.0*up + cp, 1.0*um + cm); ap = max(ap, 0.);
//   double am = min(1.0*up - cp, 1.0*um - cm); am = min(am, 0.);
//   double qm = hm*um, qp = hp*up;
//   double deltaUEps1 = (hp-hm)>0 ? max((hp-hm), tinyParam) : min((hp-hm), -1.0*tinyParam);
//   double deltaUEps2 = (qp-qm)>0 ? max((qp-qm), tinyParam) : min((qp-qm), -1.0*tinyParam);
//   double shat1 = (2.0*(qp-qm))/((hp-hm)+deltaUEps1);
//   double shat2 = (2.0*((qp*up + G*pow(hp ,2.0)/2.)-(qm*um + G*pow(hm ,2.0)/2.)))/((qp-qm)+deltaUEps2);
//   double smax = max(shat1, shat2), smin = min(shat1, shat2);
//   if (smax>tinyParam)
//   {
//     ap = min(ap, smax);
//     am = max(am, -1.0*smax);
//   }
//   else if (smin<(-1.0*tinyParam))
//   {
//     ap = min(ap, -1.0*smin);
//     am = max(am, smin);
//   }
//   double a = max(ap, -am);
//   if (a > epsilon) {
//     double wint = (ap*hp-am*hm-(qp-qm))/(ap-am);
//     double qCorr = mmo((hp-wint)/(ap-am), (wint-hm)/(ap-am));
//     *fh = (ap*qm - am*qp)/(ap - am) + (ap*am*((hp - hm)/(ap - am)-qCorr)); // (4.5) of [1]
//     wint = (ap*qp-am*qm-((1.0*qp*up + G*pow(hp ,2.0)/2.)-(1.0*qm*um + G*pow(hm ,2.0)/2.)))/(ap-am);
//     qCorr = mmo((qp-wint)/(ap-am), (wint-qm)/(ap-am));
//     *fq = (ap*(1.0*qm*um + G*pow(hp ,2.0)/2.) - am*(1.0*qp*up + G*pow(hp ,2.0)/2.))/(ap - am) + (ap*am*((qp - qm)/(ap - am)-qCorr));
//
//     //*fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
//     //*fq = (ap*(1.0*qm*um + G*sq(hm)/2.) - am*(1.0*qp*up + G*sq(hp)/2.) +
// 	    //ap*am*(qp - qm))/(ap - am);
//     double dt = CFL*Delta/a;
//     if (dt < *dtmax)
//       *dtmax = dt;
//   }
//   else
//     *fh = *fq = 0.;
// }
