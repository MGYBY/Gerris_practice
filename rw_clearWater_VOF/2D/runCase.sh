# python3 create_mesh.py
python3 create_mesh.py

# merge parameter file
cat meshFile >> column.gfs

echo "Mesh generation completed."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D -m shock_adapt1.gfs
echo "Start to run the main simulation program ... ..."
# mpirun -np 4 gerris2D -m shock_adapt1.gfs
gerris2D -m column.gfs

# gerris2D -m shock_adapt1.gfs
