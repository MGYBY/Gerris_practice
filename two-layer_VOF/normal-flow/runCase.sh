# gauge location file (numerical wave gauge, center of the domain)
awk 'BEGIN{ for (y = 0; y <= 3.3333333; y += 3.3333333/100.) print 0.00003, y, 0.;}' > centerlineGaugeLoc

echo "Gauges finished."

# python3 create_mesh.py
python3 create_mesh_parallel_safe.py

# merge parameter file
cat meshFile >> rw.gfs

echo "Mesh finished."

touch time-stats.txt
touch proj-stats.txt
touch diff-stats.txt

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main simulation program ... ..."
# gerris2D -s -m 1 rw.gfs > splitsimulation.gfs
# gerris2D -b -m 3 rw.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D -m parallelsimulation.gfs
gerris2D -m rw.gfs
# mpirun -np 4 gerris2D -m rw_adaptYS.gfs
# mpirun -np 4 gerris2D -m rw.gfs
