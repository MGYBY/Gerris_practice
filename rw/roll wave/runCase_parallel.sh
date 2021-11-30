python3 create_mesh.py

# merge parameter file
cat meshFile >> shock_adapt1.gfs

# run Gerris in parallel
gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
mpirun -np 3 gerris2D -m parallelsimulation.gfs

# gerris2D -m shock_adapt1.gfs