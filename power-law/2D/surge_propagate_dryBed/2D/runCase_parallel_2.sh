gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt

cat boxes.txt >> pl_ext_fr.gfs

echo "Mesh finished."

# gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt

# cat boxes.txt >> pl_ext_fr.gfs

touch time-stats.txt

awk 'BEGIN{
nd = 0.38093709;
nv = 3.3915176;
sl = 0.06;
grav = 9.81;
xb = 15.0/sl*nd;
width = 75.0*nd;
lbox = width*5.0;
numBoxX = 3.0
ycenter = width*5.0/2.0;
lx = numBoxX*lbox;
maxLevel = 9;
numPoint = numBoxX*(2^maxLevel);
for (x = 1; x <= numPoint; x += 1) {printf("%9f %9f %3f\n", xCoord, ycenter, 0.0); xCoord+=((lx-0.0)/(numPoint-1)); };
}' > gaugeLoc

# run Gerris in parallel
# gerris2D -m -b 3 pl_ext_fr.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs
# mpirun -np 3 gerris2D -m pl_ext_fr.gfs

# gerris2D -s 1 shock_adapt1.gfs > splitsimulation.gfs
# gerris2D -b 3 splitsimulation.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

gerris2D -m pl_ext_fr.gfs
