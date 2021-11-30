python3 create_mesh.py

# merge parameter file
cat meshFile >> shock_adapt1.gfs

# gauge location file (centreline)
awk 'BEGIN{ for (x = 0.0; x <= 42.0; x += 42.0/5000.) print x, 1.50, 0.;}' > gaugeLoc

# run Gerris
gerris2D -m shock_adapt1.gfs