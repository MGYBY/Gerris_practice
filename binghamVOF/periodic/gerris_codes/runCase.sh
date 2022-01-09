# gauge location file (numerical wave gauge, center of the domain)
awk 'BEGIN{ for (y = 0.04; y <= 0.50; y += 0.46/200.) print 9.131, y, 0.;}' > gaugeLoc

echo "Gauges finished."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main simulation program ... ..."
# gerris2D -s -m 1 rw.gfs > splitsimulation.gfs
# gerris2D -b -m 3 rw.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D -m parallelsimulation.gfs
gerris2D -m rw.gfs
# mpirun -np 4 gerris2D -m rw.gfs
