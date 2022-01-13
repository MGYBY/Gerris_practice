import numpy as np
import csv
import os

filename = "probe1"
a = np.loadtxt(filename, dtype=np.float32)

depth = a[:,1]
time = a[:,0]

#depth_filtered = depth
count = 0
depth_filtered = np.zeros_like(depth)
#depth_filtered[0] = depth[0]
#depth_filtered[1] = depth[1]

num_gauge = np.size(depth)
jump_cond = 5.0

i = 0

# traverse the data file
while i<(num_gauge-2):
#for i in range(2,num_gauge-2):
    # j index stores the length of the same time moment
    j = i
    while(time[i]==time[j]):
        # deal with OOB
        if ((j+1)<=(num_gauge-1)):
            j = j+1
        else:
            j = num_gauge-1
            break

    if (j == num_gauge-1):
        break

    seg_max = np.max(depth[i:j+1])
    seg_min = np.min(depth[i:j+1])

    if (count<=1):
        depth_filtered[count] = seg_max
    else:
    # compare the extreme value of the data segment
        if (np.abs(seg_max - depth_filtered[count]) > jump_cond*np.abs((depth_filtered[count-1] - depth_filtered[count-2]))):
            depth[i:j+1] = seg_min
        elif (np.abs(seg_min - depth_filtered[count-1]) > jump_cond*np.abs((depth_filtered[count-1] - depth_filtered[count-2]))):
            depth[i:j+1] = seg_max

    i = j+1
    count = count+1

    print("i=%d" % i)

    #if (((i)<num_gauge/4.0) and ((i+1)>num_gauge/4.0)):
        #print("25 percent finished.")

    #if (((i)<num_gauge/2.0) and ((i+1)>num_gauge/2.0)):
        #print("50 percent finished.")

    #if (((i)<3.0*num_gauge/4.0) and ((i+1)>3.0*num_gauge/4.0)):
        #print("75 percent finished.")

# write the filtered data file
file_name = "probe1_filtered"
if (not os.path.exists('./%s' % file_name)):
    with open(file_name,'w') as file:
        print("Writing filtered data to file ... ...")
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(zip(np.transpose(time),np.transpose(depth)))
        # file.write("{0} {1}\n".format(sim_time, bow_shock_standoff))
