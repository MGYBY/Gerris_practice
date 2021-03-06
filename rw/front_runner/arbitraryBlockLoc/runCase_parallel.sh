# python3 create_mesh.py
python3 create_mesh_parallel.py

# merge parameter file
cat meshFile >> shock_adapt1.gfs

echo "Mesh finished."

# gauge location file (centreline)
awk 'BEGIN{ for (x = 0.0; x <= 162.0; x += 162.0/24000.) print x, 1.500, 0.;}' > gaugeLoc

echo "Slice gauges finished."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main simulation program ... ..."
mpirun -np 4 gerris2D -m shock_adapt1.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -m shock_adapt1.gfs