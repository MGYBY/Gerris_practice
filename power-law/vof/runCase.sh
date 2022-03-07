gerris3D -s 1 -m  impact.gfs > splitsimulation.gfs
gerris3D -b 4 splitsimulation.gfs > parallelsimulation.gfs
mpirun -np 4 gerris3D -m parallelsimulation.gfs
