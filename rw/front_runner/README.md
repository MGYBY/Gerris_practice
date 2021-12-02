# The roll wave impact simulation of the front runner
1. Use generalized MINMOD limiter (need to modify source code): kurganov (r=1.30), yu (r=1.3667). Severe oscillation observed by Sweby limiter.
2. Careful selection of the adaptivity criteria: slightly strict at the beginning.
3. Output frquency control.
4. CFL number should be sufficiently small (<=0.30) to control oscillation.
5. Kinetic scheme seems to be better to control oscillation.
6. Use `GfsAdaptError` and `GfsAdaptGradient` to control oscillation.
7.  - [x] centreline slice.
