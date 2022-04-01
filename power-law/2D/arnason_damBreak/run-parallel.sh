#!/bin/sh
gerris2D -m -s 1  pl_ext.gfs > splitsimulation.gfs
gerris2D -b 4 splitsimulation.gfs > parallelsimulation.gfs
mpirun -np 4 gerris2D parallelsimulation.gfs

