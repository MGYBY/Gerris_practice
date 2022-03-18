#!/bin/sh
gerris2D -m -s 1  pl_ext.gfs > splitsimulation.gfs
gerris2D -b 2 splitsimulation.gfs > parallelsimulation.gfs
mpirun -np 2 gerris2D parallelsimulation.gfs

