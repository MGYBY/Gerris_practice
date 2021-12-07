import numpy as np

l_y = 3.0
n_y = 32
dy = l_y/n_y

filename_h = "initH"

output_h = "initH2D"

# "x h u" file
depth_1d = np.loadtxt(filename_h)
num_rows, num_cols = depth_1d.shape
with open(output_h,'w') as file:
    for j in range(0, num_rows):
        for i in range(0, n_y):
            file.write("{0} {1} {2}".format(depth_1d[j,0], ((i-0.50)*dy), depth_1d[j,1]))
            file.write('\n')


# [depth_1d[j,0], ((i-0.50)*dy), depth_1d[j,1]]