# gauge location file (numerical wave gauge, center of the domain)
# awk 'BEGIN{ for (y = 0.0; y <= 0.0042; y += 0.0042/100.) print 0.001045, y, 0.;}' > gaugeLocHeight

awk 'BEGIN {
    numBox = 3.0;
    domainLength = 0.1158174;
    normalDepth = 0.00160556;
    nx = 3.0*(2^9);
    xEnd = ((nx-1.0)/nx)*domainLength*(numBox-1.0)/numBox;
    yEnd = 2.0*normalDepth;
    numPoint = 120.0;
    delta = yEnd/numPoint;
    for (y = yEnd/numPoint/2.0; y <= yEnd; y += yEnd/numPoint) print xEnd, y, 0.;
}' > gaugeLoc1

awk 'BEGIN {
    numBox = 3.0;
    domainLength = 0.1158174;
    normalDepth = 0.00160556;
    nx = 3.0*(2^9);
    xEnd = domainLength*(1.0/nx+(-1.0)/numBox);
    yEnd = 2.0*normalDepth;
    numPoint = 120.0;
    delta = yEnd/numPoint;
    for (y = yEnd/numPoint/2.0; y <= yEnd; y += yEnd/numPoint) print xEnd, y, 0.;
}' > gaugeLoc2

echo "Gauges finished."

# touch totalEnergy
# touch kineticEnergy
# touch potentialEnergy
touch ef # energy flux file

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main simulation program ... ..."
gerris2D -m -s 3 rw.gfs > splitsimulation.gfs
gerris2D -m -b 4 splitsimulation.gfs > parallelsimulation.gfs
mpirun -np 4 gerris2D -m parallelsimulation.gfs
# gerris2D -m rw.gfs

# gerris2D -m -b 4 rw.gfs > parallelsimulation.gfs
# mpirun -np 4 gerris2D -m parallelsimulation.gfs

# mpirun -np 4 gerris2D -m rw_adaptYS.gfs
# mpirun -np 6 gerris2D -m rw.gfs
