Modified Gerris codes to add the functionality.

1. Added "Kurganov" feature: avoid the complication of the new eigenstructure by using central-upwind scheme. For SWE and SW PL fluid.
2. Added "slopeBoolean" feature: distinguish X- and Y- sweep of the dimension-by-dimension approach. <0.50: slope coordinate; >0.50: normal SWE formulation.
3. Added "yu_limiter" feature: tunable generalized MINMOD limiter.
4. Timestep determination based on eigenvalues of power-law system.
5. Added "MC" limiter. Improved resolution.
6. Added SGVL and SGVA limiters.
7. Added van-Leer with Klein's sharpening.
8. Added central-upwind scheme with numerical dissipation reduction (following Kurganov & Lin (2007)).
