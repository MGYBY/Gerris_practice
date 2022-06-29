# a more reasonable way to automate the guage generation
awk 'BEGIN {
    numBox = 5.0;
    normalVel = 0.0456607;
    normalDepth = 0.00151404;
    gravCoeff = 9.81;
    channelSin = 0.06;
    domainLength = 50.4903*(normalVel^2)/(gravCoeff*channelSin);
    nx = 5.0*(2^9);
    xEnd = ((nx-1.0)/nx)*domainLength*(numBox-1.0)/numBox;
    yEnd = 2.0*normalDepth;
    numPoint = 80.0;
    delta = yEnd/numPoint;
    for (y = yEnd/numPoint/2.0; y <= yEnd; y += yEnd/numPoint) print xEnd, y, 0.;
}' > gaugeLoc1

awk 'BEGIN {
    numBox = 5.0;
    normalVel = 0.0456607;
    normalDepth = 0.00151404;
    gravCoeff = 9.81;
    channelSin = 0.06;
    domainLength = 50.4903*(normalVel^2)/(gravCoeff*channelSin);
    nx = 5.0*(2^9);
    xEnd = domainLength*(1.0/nx+(-1.0)/numBox);
    yEnd = 2.0*normalDepth;
    numPoint = 80.0;
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
