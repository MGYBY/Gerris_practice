#!/bin/bash
# python3 create_mesh.py
python3 create_mesh_parallel.py

# merge parameter file
cat meshFile >> shock_adapt1.gfs

# gauge location file (centreline)
awk 'BEGIN{ for (x = 0.0; x <= 162.0; x += 162.0/24000.) print x, 1.500, 0.;}' > gaugeLoc

awk 'BEGIN{
cellSize = 3.0/512.0
spacing = cellSize/4.0
leftXCoord = 160.0-spacing*2.0
leftYCoord = 1.50
rightXCoord = 160.0+0.30*(3^0.50)
rightYCoord = 1.80+spacing/(3^0.50)*2.0
numPoint = 120
xCoord = leftXCoord
yCoord = leftYCoord
for (x = 1 ; x <= (numPoint); x += 1) {printf("%9f %9f %3f\n", xCoord, yCoord, 0.0); xCoord+=((rightXCoord-leftXCoord)/(numPoint-1)) ; yCoord+=((rightYCoord-leftYCoord)/(numPoint-1)); }
}' > RGaugeLoc

# run Gerris in parallel
# gerris2D -m -b 3 shock_adapt1.gfs > parallelsimulation.gfs
mpirun -np 5 gerris2D -m shock_adapt1.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -m shock_adapt1.gfs
