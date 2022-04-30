# The attempt to reproduce Liu & Mei 1994 results using NS-VOF method.

Useful links:

1. http://basilisk.fr/sandbox/vatsal/GenaralizedNewtonian/SurfaceWavesBingham.c
2. http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_simple.c
3. http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_cohesif.c
4. http://basilisk.fr/sandbox/vatsal/GenaralizedNewtonian/Couette_NonNewtonian.c
5. http://basilisk.fr/sandbox/lucaspagoto/Elasto_viscoplastic/shear_droplet2.c
6. http://basilisk.fr/sandbox/M1EMN/Exemples/column_SCC.c
7. http://basilisk.fr/sandbox/hiranya/viscoplasticsheet.c
8. http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_tube.c
9. Some free-surface treatments:

        http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c

11. Refer to these codes for inclusion of gravity: 

                                                http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_cohesif.c

                                                http://basilisk.fr/src/test/poiseuille45.c

11. Solid force: 

                 http://basilisk.fr/sandbox/ghigo/src/test-navier-stokes/cylinder-accelerating.c

                 http://basilisk.fr/sandbox/ghigo/src/test-navier-stokes/cylinder-strouhal.c
                 
                 http://basilisk.fr/sandbox/ghigo/src/myspheroid.h
                 
                 http://basilisk.fr/src/test/starting.c
                 
                 http://gerris.dalembert.upmc.fr/gerris/examples/examples/ship.html#htoc17
                 
12. Gerris examples:

```
http://gfs.sourceforge.net/wiki/index.php/Cremonesi_et_al_benchmark
http://gerris.dalembert.upmc.fr/gerris/tests/tests/couette.html
http://gerris.dalembert.upmc.fr/gerris/examples/examples/column.html#htoc9
```

13. Gerris VTK or Tecplot format output for parallel running:

```
https://sourceforge.net/p/gfs/mailman/message/27056583/
```

14. Trick to output VOF surface every certain time interval:
```
http://gfs.sourceforge.net/wiki/index.php/Gfsview-batch
```

15. "Comet droplet" treatment:
```
https://groups.google.com/g/basilisk-fr/c/HO_u2ErNhO0/m/3w42XrajAAAJ
```
