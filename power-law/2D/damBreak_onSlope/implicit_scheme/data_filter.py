import numpy as np
import csv
import os

filename = "frontPos"
a = np.loadtxt(filename, dtype=np.float32)

x_coord = a[:,1]
time = a[:,0]

#depth_filtered = depth
count = 0
x_filtered = np.zeros_like(x_coord)
time_filtered = np.zeros_like(time)
#depth_filtered[0] = depth[0]
#depth_filtered[1] = depth[1]

num_gauge = np.size(x_coord)

i = 0

# traverse the data file
while i<(num_gauge-2):
#for i in range(2,num_gauge-2):
    # j index stores the length of the same time moment
    j = i+1

    if (count<=0):
        x_filtered[count] = x_coord[i]
        time_filtered[count] = time[i]
        count += 1
    else:
    # compare the extreme value of the data segment
        while(j<=num_gauge-1):
            if (np.abs(x_coord[i] - x_coord[j]) > 1.e-8):
                x_filtered[count] = x_coord[i]
                time_filtered[count] = time[i]
                # j += 1
                count += 1
                break
            else:
                j += 1
                i += 1
        i += 1

time_filtered_mod = time_filtered[np.where(time_filtered>0)]
x_filtered_mod = x_filtered[np.where(x_filtered>0)]

# write the filtered data file
file_name = "frontPos_filtered"
if (not os.path.exists('./%s' % file_name)):
    with open(file_name,'w') as file:
        print("Writing filtered data to file ... ...")
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(zip(np.transpose(time_filtered_mod),np.transpose(x_filtered_mod)))
        # file.write("{0} {1}\n".format(sim_time, bow_shock_standoff))
