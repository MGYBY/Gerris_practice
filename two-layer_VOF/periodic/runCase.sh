# for energy contents
# touch totalEnergy kineticEnergy potentialEnergy

gerris2D -s 3 -m  rw.gfs > splitsimulation.gfs
gerris2D -b 9 splitsimulation.gfs > parallelsimulation.gfs
# # mpirun -np 4 gerris3D -m parallelsimulation.gfs
mpirun --oversubscribe -np 9 gerris2D parallelsimulation.gfs -parallel

# gerris2D -m rw.gfs
