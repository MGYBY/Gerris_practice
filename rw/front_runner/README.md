# The roll wave impact simulation of the front runner
1. Use generalized MINMOD limiter (need to modify source code): kurganov (r=1.30), yu (r=1.3667). Severe oscillation observed by Sweby limiter. Another cure: lower MINLEVEL to 3.
2. Careful selection of the adaptivity criteria: slightly strict at the beginning.
3. Output frquency control.
4. CFL number should be sufficiently small (<=0.30) to control oscillation.
5. Kinetic scheme seems to be better to control oscillation.
6. Use `GfsAdaptError` and `GfsAdaptGradient` to control oscillation.
7.  - [x] centreline slice.
8. Sweby limiter has **large time-step restriction** compared to minmod at (some) high Fr (for example, Fr=5.60). Reason unknown. Possible reason: unmatched limiter functions for embedded solid and internal fluid field. A better option is to use van-leer limiter:
```
gradient = gfs_center_van_leer_gradient
```
