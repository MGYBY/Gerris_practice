# python3 create_mesh.py
python3 create_mesh_parallel.py

# merge parameter file
cat meshFile >> shock_adapt1.gfs

echo "Mesh finished."

# gauge location file (centreline)
awk 'BEGIN{ for (x = 0.0; x <= 162.0; x += 162.0/24000.) print x, 0.75, 0.;}' > gaugeLoc

echo "Slice gauges finished."

# mapped initial input file
python3 format_text.py
sh init.sh

echo "Mapped initial condition finished."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D -m shock_adapt1.gfs
echo "Beginning main simulation program... ..."
mpirun -np 4 gerris2D -m shock_adapt1.gfs

# gerris2D -m shock_adapt1.gfs
