# gauge location file (numerical wave gauge, center of the domain)
awk 'BEGIN{ for (y = 0.0013; y <= 0.0083; y += 0.007/100.) print 0.0, y, 0.;}' > gaugeLoc1
awk 'BEGIN{ for (y = 0.0013; y <= 0.0083; y += 0.007/100.) print 0.0, y, 0.;}' > gaugeLoc2
awk 'BEGIN{ for (y = 0.0013; y <= 0.0083; y += 0.007/100.) print 0.0, y, 0.;}' > gaugeLoc3
awk 'BEGIN{ for (y = 0.0013; y <= 0.0083; y += 0.007/100.) print 0.0, y, 0.;}' > gaugeLoc4
awk 'BEGIN{ for (y = 0.0013; y <= 0.0083; y += 0.007/100.) print 0.0, y, 0.;}' > gaugeLoc5
echo "Gauges finished."

# python3 create_mesh.py
python3 create_mesh_parallel_safe.py
# merge parameter file
cat meshFile >> rw.gfs
echo "Mesh finished."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main simulation program ... ..."
# gerris2D -s -m 1 rw.gfs > splitsimulation.gfs
# gerris2D -b -m 3 rw.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D -m parallelsimulation.gfs
# gerris2D -m rw.gfs
# mpirun -np 4 gerris2D -m rw_adaptYS.gfs
mpirun -np 6 gerris2D -m rw.gfs
