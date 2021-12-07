import numpy as np

l_y = 3.0
n_y = 32
dy = l_y/n_y

filename_h = "initH"

output_h = "initH2D"
output_u = "initU2D"

# "x h u" file
depth_1d = np.loadtxt(filename_h)
num_rows, num_cols = depth_1d.shape
with open(output_h,'w') as file:
    for j in range(0, num_rows):
        # first take care of the upper boundary
        file.write("{0} {1} {2}".format(depth_1d[j,0], 0.0, depth_1d[j,1]))
        file.write('\n')
        for i in range(0, n_y):
            file.write("{0} {1} {2}".format(depth_1d[j,0], ((i+1.0-0.50)*dy), depth_1d[j,1]))
            file.write('\n')
        # lastly, take care of the lower boundary
        file.write("{0} {1} {2}".format(depth_1d[j,0], l_y, depth_1d[j,1]))
        file.write('\n')

with open(output_u,'w') as file:
    for j in range(0, num_rows):
        # first take care of the upper boundary
        file.write("{0} {1} {2}".format(depth_1d[j,0], 0.0, depth_1d[j,2]))
        file.write('\n')
        for i in range(0, n_y):
            file.write("{0} {1} {2}".format(depth_1d[j,0], ((i+1.0-0.50)*dy), depth_1d[j,2]))
            file.write('\n')
        # lastly, take care of the lower boundary
        file.write("{0} {1} {2}".format(depth_1d[j,0], l_y, depth_1d[j,2]))
        file.write('\n')


# [depth_1d[j,0], ((i-0.50)*dy), depth_1d[j,1]]