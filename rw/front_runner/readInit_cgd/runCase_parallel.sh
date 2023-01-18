# python3 create_mesh.py
python3 create_mesh_parallel.py

# merge parameter file
cat meshFile >> shock_adapt1.gfs

echo "Mesh finished."

./constructCGD.sh

echo "CGD files finished."

# gauge location file (centreline)
awk 'BEGIN {
    nd = 0.00798;
    numBox = 1004.0;
    maxLevel = 3;
    lx = 40.596;
    ly = nd*1.0*5.0;
    lCenter = ly/2.0;
    numPoint = numBox*(2^maxLevel);
    for (x = lx/numPoint/2.0; x <= lx; x += lx/numPoint) print x, lCenter, 0.;
}' > gaugeSlice

awk 'BEGIN {
    channelSin = 0.05011*1.0;
    froude = 3.71;
    normalDepth = 0.00798;
    gravCoeff = 9.81;
    ly = normalDepth*1.0*5.0;
    normalVel = froude*((gravCoeff*((1.0-channelSin*channelSin)^0.50)*normalDepth)^0.50);
    maxLevel = 3;
    numBox = 1004.0;
    numCell = numBox*(2^maxLevel);
    xCoord = 251.1779*normalDepth/channelSin*((numCell-1.00)/numCell);
    yCenter = ly/2.0;
    printf "%g %g 0.", xCoord, yCenter;
}' > gaugeRmax

# exit

echo "Slice gauges finished."

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
echo "Beginning main program ... ..."
mpirun -np 6 gerris2D -m shock_adapt1.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -m shock_adapt1.gfs
