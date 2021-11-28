python3 create_mesh.py

# merge parameter file
cat meshFile >> shock.gfs

# run Gerris
gerris2D -m shock.gfs