1D dam-break on slope, power-law fluids.

Features:
- Central-upwind scheme based shock-capturing scheme. Dimensional form for now.
- Backward linearized Euler scheme for the source term.
- Bi-tree adaptivity based on **(1) depth gradient (2) depth and (3) Froude-number gradient**.
- Wave-front tracking (adopt the more efficient way by M1EMN), extent-averaged depth and velocity. Filtered wavefront (using `data_filter.py`).
- Global maximum depth and velocity. Maximum $\rm Fr$.
- Text output and image output. (TODO: output $\rm Fr$)
- Initial-stage timestep ($\Delta t_{max}$) restriction to better capture the wavefront and celerity.
