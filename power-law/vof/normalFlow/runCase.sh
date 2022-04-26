# gauge location file (numerical wave gauge, center of the domain)
awk 'BEGIN{ for (y = 0; y <= 0.002791; y += 0.002791/32.) print 0.0020933486, y, 0.;}' > centerlineGaugeLoc

echo "Gauges finished."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main simulation program ... ..."
# gerris2D -s -m 1 rw.gfs > splitsimulation.gfs
# gerris2D -b -m 3 rw.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D -m parallelsimulation.gfs
gerris2D -m rw.gfs
# mpirun -np 4 gerris2D -m rw_adaptYS.gfs
# mpirun -np 4 gerris2D -m rw.gfs
