gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt

cat boxes.txt >> pl_ext_fr.gfs

touch time-stats.txt

awk 'BEGIN{
normalDepth = 0.002099551886;
normalVel = 0.1433858;
grav = 9.81;
slope = 0.060;
xb = 50.0*normalVel^2/(slope*grav);
width = normalDepth*50.0;
yCenterCoord = width*5.0/2.0;
cellSize = width*5.0/(3.0*128.0);
spacing = cellSize/4.0;
leftXCoord = xb-spacing*1.75;
leftYCoord = yCenterCoord;
rightXCoord = xb+1.0*width;
rightYCoord = yCenterCoord+(width/2.0)+spacing/(3^0.50)*1.75;
numPoint = 120;
xCoord = leftXCoord;
yCoord = leftYCoord;
for (x = 1 ; x <= (numPoint); x += 1) {printf("%9f %9f %3f\n", xCoord, yCoord, 0.0); xCoord+=((rightXCoord-leftXCoord)/(numPoint-1)) ; yCoord+=((rightYCoord-leftYCoord)/(numPoint-1)); }
}' > RGaugeLoc

# run Gerris in parallel
gerris2D -m -b 3 pl_ext_fr.gfs > parallelsimulation.gfs
mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -m shock_adapt1.gfs
