1D dam-break on slope, power-law fluids.

Features:
- Explicit scheme: central-upwind scheme based shock-capturing scheme. Dimensional form for now.
- Implicit scheme: all-mach method.
- Backward linearized Euler scheme for the source term. Useful formula:

$\Large u_{N+1}=\displaystyle \frac{u_N+\Delta t S_o g'}{1+\Delta t \frac{\mu_n}{\rho}(\frac{1+2n}{n})^n \frac{u_N^{n-1}}{h_N^{n+1}}}$

$\large (h,\,q)$ formulation useful for implicit Saint-Venant solver and Gerris solver

$\Large q_{N+1}=\displaystyle \frac{q_N+\Delta t g' S_o h_N}{1+\Delta t \frac{\mu_n}{\rho} (\frac{1+2n}{n})^n \frac{q_N^{n-1}}{h_N^{2 n}}}$
- Bi-tree adaptivity based on **(1) depth gradient (2) depth and (3) Froude-number gradient**.
- Wave-front tracking (adopt the more efficient way by M1EMN), extent-averaged depth and velocity. Filtered wavefront (using `data_filter.py`).
- Global maximum depth and velocity. Maximum $\rm Fr$.
- Energy content: potential energy and kinetic energy.
- Text output and image output. $x$ | $h$ | $u$ | $\rm Re$ | $\rm Fr$
- Initial-stage timestep ($\Delta t_{max}$) restriction to better capture the wavefront and celerity.
