# gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt
#
# cat boxes.txt >> pl_ext_fr.gfs

touch time-stats.txt

awk 'BEGIN{
level = 9;
numBoxX = 6;
numPoint = numBoxX*(2^level);
initDepth = 4.0;
aspectRatio = 5.0;
obsWidth = aspectRatio*initDepth;
boxLength = 5.0*obsWidth;
ycenter = boxLength/2.0;
lx = boxLength*numBoxX;
xb = 12.0*initDepth*aspectRatio;
xCoord = xb - lx/numPoint/2.0;
printf("%9f %9f %3f\n", xCoord, ycenter, 0.0);
# for (x = 1; x <= numPoint; x += 1) {printf("%9f %9f %3f\n", xCoord, ycenter, 0.0); xCoord+=((lx-0.0)/(numPoint-1)); };
}' > runupLoc

awk 'BEGIN{
level = 9;
numBoxX = 6;
numPoint = numBoxX*(2^level);
initDepth = 4.0;
aspectRatio = 5.0;
obsWidth = aspectRatio*initDepth;
boxLength = 5.0*obsWidth;
ycenter = boxLength/2.0;
lx = boxLength*numBoxX;
xb = 12.0*initDepth*aspectRatio;
xCoord = 0.0;
for (x = 1; x <= numPoint; x += 1) {printf("%9f %9f %3f\n", xCoord, ycenter, 0.0); xCoord+=((lx-0.0)/(numPoint-1)); };
}' > gaugeLoc

# run Gerris in parallel
# gerris2D -m -b 3 pl_ext_fr.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs
# mpirun -np 3 gerris2D -m pl_ext_fr.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

gerris2D -m pl_ext_db.gfs
