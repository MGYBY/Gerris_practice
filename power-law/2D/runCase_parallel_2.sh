gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt

cat boxes.txt >> pl_ext_fr.gfs

touch time-stats.txt

# run Gerris in parallel
gerris2D -m -b 3 pl_ext_fr.gfs > parallelsimulation.gfs
mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -m shock_adapt1.gfs
