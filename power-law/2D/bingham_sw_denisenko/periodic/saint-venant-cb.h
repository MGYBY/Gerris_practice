/**
# A solver for the Two-layer equations
The Saint-Venant equations for two layers:

The [Saint-Venant equations](http://en.wikipedia.org/wiki/Shallow_water_equations)
can be written in integral form as the hyperbolic system of
conservation laws
$$
\partial_t \int_{\Omega} \mathbf{q} d \Omega +
\int_{\partial \Omega} \mathbf{f} (
\mathbf{q}) \cdot \mathbf{n}d \partial
\Omega + \int_{\Omega} \mathbf{F} d \Omega= 0
$$
where $\Omega$ is a given subset of space, $\partial \Omega$ its boundary and
$\mathbf{n}$ the unit normal vector on this boundary. For
conservation of mass and momentum in the shallow-water context, $\Omega$ is a
subset of bidimensional space and $\mathbf{q}$ and
$\mathbf{f}$ are written
$$
\mathbf{q} =
h_1,	\;
h_1u_1_x,\;
h_1u_1_y,\;
h_2 ,\;
h_2u_2_x,\;
h_2u_2_y
$$
\;\;\;\;\;\;
$$
\mathbf{f} (\mathbf{q}) = \left(\begin{array}{cc}
h_1u_1_x & h_1u_1_y\					\
h_1u_1_x^2 + \frac{1}{2} gh_1^2 & h_1u_1_xu_1_y\	\
h_1u_1_xu_1_y & h_1u_1_y^2 + \frac{1}{2} gh_1^2\	\
h_2u_2_x & h_2u_2_y\					\
h_2u_2_x^2 + \frac{1}{2} gh_2^2 & h_2u_2_xu_2_y\	\
h_2u_2_xu_2_y & h_2u_2_y^2 + \frac{1}{2} gh_2^2
\end{array}\right),
\;\;\;\;\;\;
\mathbf{F} = \left(\begin{array}{c}
0\								\
g h_1 \frac{\partial}{\partial x} \left( h_2 + z_b \right)\	\
g h_1 \frac{\partial}{\partial y} \left( h_2 + z_b \right)\	\
0\									\
g h_2 \frac{\partial}{\partial x} \left( z_b + \frac{\rho_1}{\rho_2} h_1 \right)\ \
g h_2 \frac{\partial}{\partial y} \left( z_b + \frac{\rho_1}{\rho_2} h_1 \right)
\end{array} \right),
$$
where $h_1$ the water depth, $h_2$ is the landslide thickness, $\mathbf{u}_1$ is the velocity vector of the water, $\mathbf{u}_2$ is the velocity vector of the landslide, and
$z_b$ the height of the topography. See also [Popinet,
2011](/src/references.bib#popinet2011) for a more detailed
introduction.

## User variables and parameters

The primary fields are the water depth $h_1$, the landslide thickness $h_2$, the bathymetry $z_b$ and
the flow speeds of water and landslide respectively $\mathbf{u}_1$ and $\mathbf{u}_2$. $\eta$ is the water level i.e. $z_b +
h_1 + h_2$. Note that the order of the declarations is important as $z_b$
needs to be refined before $h_2$, before $h_1$, before $\eta$. */

// TODO: change these weird names orginated from the two-layer model
// TODO: many function names are very outdated, make them compatible with the latest version of Basilisk

scalar zb[], hie[], eta[], h[];
vector u[];
//double default_sea_level=0.;

/**
The only physical parameter is the acceleration of gravity `G` and the densities of the two fluids.
Cells are considered "dry" when the water depth is less than the `dry` parameter (this
should not require tweaking). */

double G = 1.;
double dry = 1e-10;

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
The generic time-integration scheme in predictor-corrector.h needs
to know which fields are updated. */

scalar * evolving = {h, u, hie};

/**
We need to overload the default *advance* function of the
predictor-corrector scheme, because the evolving variables ($h_1$, $\mathbf{u}_1$, $h_2$ and
$\mathbf{u}_2$) are not the conserved variables $h_1$, $h_1\mathbf{u}_1$, $h_2$ and
$h_2\mathbf{u}_2$. */

static void advance_two_layer (scalar * output, scalar * input,
			       scalar * updates, double dt)
{
  //  fprintf(stderr,"advance_two_layer\n");
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0], dh = updates[0];
  scalar hiei = input[1], hieo = output[1], dhie = updates[1];
  vector ui = { input[2], input[3] }, uo = { output[2], output[3] };
  vector dhu = { updates[2], updates[3] };
  // new fields in ho[], uo[]
  foreach() {
    double hold = hi[];
    double hieold = hiei[];
    ho[] = hold + dt*dh[];
    hieo[] = hieold + dt*dhie[];
    eta[] = ho[] + zb[];

    if (ho[] > dry)
    {
      //      fprintf(stderr," H1 wet\n");
      foreach_dimension()
        uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];
    }
    else
    {
      //      fprintf(stderr," H1 dry\n");
      foreach_dimension()
	uo.x[] = 0.;
    }
  }
  //  fprintf(stderr,"boundary\n");
  //h1o.prolongation = refine_linear;
  boundary ({ho, eta, uo, hieo});
  //  fprintf(stderr,"done\n");
}

/**
When using an adaptive discretisation (i.e. a quadtree)., we need
to make sure that $\eta$ is maintained as $z_b + h_1 + h_2$ whenever cells are
refined or coarsened. */

/**
#if QUADTREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h1[] + h2[];
}

static void coarsen_eta (Point point, scalar eta)
{
  eta[] = zb[] + h1[] + h2[];
}
#endif
*/

/**
### Computing fluxes

Various approximate Riemann solvers are defined in [riemann.h](). */

#include "./riemann-cb.h"

double update_two_layer (scalar * evolving, scalar * updates, double dtmax)
{

/**
We first recover the currently evolving fields (as set by the
predictor-corrector scheme). */

  scalar h = evolving[0];
  scalar hie = evolving[1];
  vector u = { evolving[2], evolving[3] };

/**
`Fh1`, `Fh2` `Fq1` and `Fq2` will contain the fluxes for $h_1$, $h_2$, $h_1\mathbf{u}_1$ and $h_2\mathbf{u}_2$
respectively and `S1` and `S2` are necessary to store the asymmetric topographic
source terms. */

  vector Fh[], Fhie[], S[];
  tensor Fq[];

/**
The gradients are stored in locally-allocated fields. First-order
reconstruction is used for the gradient fields. */

  vector gh[], geta[], ghie[];
  tensor gu[];
  for (scalar s in {gh, geta, ghie, gu}){
    s.gradient = zero;
#if QUADTREE   //Is this where the problem is?  Do I need to refine some other way?
      s.prolongation = refine_linear;
    #endif
  }
  gradients ({h, eta, hie, u}, {gh, geta, ghie, gu});
  //fprintf(stdout,"%% updates\n");
/**
The faces which are "wet" on at least one side are traversed.
First we see whether "wet" in bottom fluid, if so look for lake at rest solution
$h_2+z_b=C_2$ $h_1=C_1$
If $h_2$ is dry look for lake at rest solution
$h_2=0$ $h_1+z_b=C_0$
*/

  foreach_face (reduction (min:dtmax)) {
    double hi = h[], hn = h[-1,0];
    if (hi > dry || hn > dry) {
      //      fprintf(stderr,"H2 wet");

/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */

      double dx = Delta/2.;
      double zi = eta[] - hi - h[];
      double zl = zi - dx*(geta.x[] - gh.x[]- gh.x[]);
      double zn = eta[-1,0] - hn - h[-1,0];
      double zr = zn + dx*(geta.x[-1,0] - gh.x[-1,0] - gh.x[-1,0]);
      double zlr = max(zl, zr);

      double hl = hi - dx*gh.x[];
      double up = u.x[] - dx*gu.x.x[];
      double hp = max(0., hl + zl - zlr);

      double hr = hn + dx*gh.x[-1,0];
      double um = u.x[-1,0] + dx*gu.x.x[-1,0];
      double hm = max(0., hr + zr - zlr);

      double hiem = hie[-1, 0] + dx*ghie.x[-1,0];
      double hiep = hie[] - dx*ghie.x[];

/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

      double fh, fu, fhie, fv;
      kurganov (hm, hp, um, up, hiem, hiep, Delta*cm[]/fm.x[], &fh, &fu, &fhie, &dtmax);
      // kurganovSharp (hm, hp, um, up, hiem, hiep, Delta*cm[]/fm.x[], &fh, &fu, &fhie, &dtmax);
      fv = (fh > 0. ? u.y[-1,0] + dx*gu.y.x[-1,0] : u.y[] - dx*gu.y.x[])*fh;

/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if QUADTREE
      if (is_prolongation(cell)) {
	hi = coarse(h,0,0,0);
	zi = coarse(zb,0,0,0);
      }
      if (is_prolongation(neighbor(-1,0))) {
	hn = coarse(h,-1,0,0);
	zn = coarse(zb,-1,0,0);
      }
#endif

      // double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
      // double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));

/**
#### Flux update */

      Fh.x[]   = fm.x[]*fh;
      Fhie.x[] = fm.x[]*fhie;
      // Fq.x.x[] = fm.x[]*(fu - sl);
      // S.x[]    = fm.x[]*(fu - sr);
      Fq.x.x[] = fm.x[]*(fu);
      S.x[]    = fm.x[]*(fu);
      Fq.y.x[] = fm.x[]*fv;
      //fprintf(stdout,"2 %g %g %g %g %g\n",x,y,fh,fu,fv);
    }
    else
      Fh.x[] = Fhie.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;
  }
  boundary_flux ({Fh, Fhie, S, Fq});

/**
#### Updates for evolving quantities

We store the divergence of the fluxes in the update fields. Note that
these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */

  scalar dh = updates[0], dhie = updates[1];
  vector dhu = { updates[2], updates[3] };

  foreach() {
    dh[] = (Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/(cm[]*Delta);
    dhie[] = (Fhie.x[] + Fhie.y[] - Fhie.x[1,0] - Fhie.y[0,1])/(cm[]*Delta);
    foreach_dimension(){
      dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta);
    }

    /**
    We also need to add the metric terms. They can be written (see
    eq. (8) of [Popinet, 2011](references.bib#popinet2011))
    $$
    S_g = h \left(\begin{array}{c}
    0\								\
    \frac{g}{2} h \partial_{\lambda} m_{\theta} + f_G u_y\	\
    \frac{g}{2} h \partial_{\theta} m_{\lambda} - f_G u_x
    \end{array}\right)
    $$
    with
    $$
    f_G = u_y \partial_{\lambda} m_{\theta} - u_x \partial_{\theta} m_{\lambda}
    $$

#### CHECK THIS - NOT SURE OF IT ####
    */

    double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
    double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
    double fG = u.y[]*dmdl - u.x[]*dmdt;
    dhu.x[] += h[]*(G*h[]/2.*dmdl + fG*u.y[]);
    dhu.y[] += h[]*(G*h[]/2.*dmdt - fG*u.x[]);

  }

  //printf ("%% Flux_updates t=%g\n",t);
  //output_field ({h1, h2, zb, Fh1x,Fh1y,Fh2x,Fh2y,S1x,S1y,S2x,S2y,dh1,dh2}, stdout, n = 1 << 10, linear = false);
  return dtmax;
}
/**
We use the main time loop (in the predictor-corrector scheme) to setup
the initial defaults. */

event defaults (i = 0)
{

/**
We overload the default 'advance' and 'update' function of the predictor-corrector
scheme and setup the refinement and coarsening methods on quadtrees. */

  advance = advance_two_layer;
  update = update_two_layer;
#if QUADTREE
  for (scalar s in {h,zb,u,hie}) {
    s.refine = s.prolongation = refine_linear;
    s.coarsen = coarsen_volume_average;
  }
  //eta.refine  =refine_eta;
  //eta.prolongation = refine_linear;
  //eta.coarsen = coarsen_eta;
#endif
}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  foreach()
    eta[] = zb[] + h[];
  boundary (all);
}


/**
## Conservation of water surface elevation

When using the default adaptive reconstruction of variables, the
Saint-Venant solver will conserve the water depth when cells are
refined or coarsened. However, this will not necessarily ensure that
the "lake-at-rest" condition (i.e. a constant water surface elevation)
is also preserved. In what follows, we redefine the `refine()` and
`coarsen()` methods of the water depth $h$ so that the water surface
elevation $\eta$ is conserved.

We start with the reconstruction of fine "wet" cells: */
#if QUADTREE
static void refine_elevation1 (Point point, scalar h1)
{
  struct { double x, y; } g; // gradient of eta
  if (h[] >= dry) {
    double eta = zb[] + h[];   // water surface elevation
    foreach_dimension()
      g.x = gradient (zb[-1,0] + h[-1,0], eta, zb[1,0] + h[1,0])/4.;
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta + g.x*child.x + g.y*child.y - zb[]);
  }
  else {
/**
The "dry" case is a bit more complicated. We look in a 3x3
neighborhood of the coarse parent cell and compute a depth-weighted
average of the "wet" surface elevation $\eta$. We need to do this
because we cannot assume a priori that the surrounding wet cells are
necessarily close to e.g. $\eta = 0$. */
    double v = 0., eta1 = 0.; // water surface elevation
    // 3x3 neighbourhood
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h1[i,j] >= dry) {
	  eta1 += h[i,j]*(zb[i,j] + h[i,j]);
	  v += h[i,j];
	}
    if (v > 0.)
      eta1 /= v; // volume-averaged eta of neighbouring wet cells
    else
      /**
If none of the surrounding cells is wet, we assume a default sealevel
at zero.
XXXXXXXXXXXX
Is this best?  What if zb is negative?  Should we set eta to that?*/
      eta1 = 0.;
/**
We then reconstruct the water depth in each child using $\eta$ (of the
parent cell i.e. a first-order interpolation in contrast to the wet
case above) and $z_b$ of the child cells. */
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta1 - zb[]);
  }
}

/**
Cell coarsening is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */
static void coarsen_elevation1 (Point point, scalar h1)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h1[] > dry) {
      eta += h[]*(zb[] + h[]);
      v += h[];
    }
/**
... and use this in combination with $z_b$ (of the coarse cell) to
compute the water depth $h$.  */
  if (v > 0.)
    h[] = max(0., eta/v - zb[]);
  else // dry cell
    h[] = 0.;
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  h.refine  = refine_elevation1;
  h.prolongation  = refine_linear;
  h.coarsen = coarsen_elevation1;
}
#else // Cartesian
void conserve_elevation (void) {}
#endif

/**
## "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (`ref`). */

#define radiation1(ref) (sqrt (G*max(h1[],0.)) - sqrt(G*max((ref) - zb[], 0.)))
#define radiation2(ref) (sqrt (G*max(h2[],0.)) - sqrt(G*max((ref) - zb[], 0.)))

/**
## Tide gauges

An array of `Gauge` structures passed to `output_gauges()` will create
a file (called `name`) for each gauge. Each time `output_gauges()` is
called a line will be appended to the file. The line contains the time
and the value of each scalar in `list` in the (wet) cell containing
`(x,y)`. The `desc` field can be filled with a longer description of
the gauge. */

  typedef struct {
    char * name;
    double x, y;
    char * desc;
    FILE * fp;
  } Gauge;

void output_gauges (Gauge * gauges, scalar * list)
{
  for (Gauge * g = gauges; g->name; g++) {
    if (!g->fp) {
      g->fp = fopen (g->name, "w");
      if (g->desc)
	fprintf (g->fp, "%s\n", g->desc);
    }
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    Point point = locate (xp, yp);
    if (point.level >= 0 && h[] > dry) {
      fprintf (g->fp, "%g", t);
      for (scalar s in list)
	fprintf (g->fp, " %g", s[]);
    }
    else {
      fprintf (g->fp, "%g", t);
      for (scalar s in list)
	fprintf (g->fp, " NaN");
    }
    fputc ('\n', g->fp);
    fflush (g->fp);
  }
}
