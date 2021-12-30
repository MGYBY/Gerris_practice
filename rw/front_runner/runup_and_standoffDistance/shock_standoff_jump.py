#!/usr/bin/env python

import numpy as np
# import time
import sys
# import re
import fnmatch
import os

# normal_depth = 0.00798
obs_loc = 160.00
sim_time = sys.argv[1]
# print('Simulation time is:{0}'.format(sim_time))
# a two-column file
# a = np.loadtxt("centerlineSlice", usecols = (1,4), skiprows=1, dtype=np.float32)
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, 'centerlineSlice*'):
        a = np.loadtxt(file, usecols = (1,4), skiprows=1, dtype=np.float32)

val = a[:,1]
loc = a[:,0]

num_gauge = np.size(val)
max_grad = 0.0
max_grad_loc = loc[0]

# locate the first discontinuity
for i in range(1,num_gauge-1):
    if((np.abs(val[i+1]-val[i])>1.99*(np.abs(val[i]-val[i-1]))) and (np.abs(obs_loc-loc[i])>0.01)):
        max_grad = val[i+1]-val[i-1]
        max_grad_loc = loc[i]
        break

bow_shock_standoff = obs_loc-max_grad_loc+(3.0/512.0)*0.55
# print('Location of bow shock is:{0}'.format(bow_shock_loc))
with open('bowShockWidth','a') as file:
    file.write("{0} {1}\n".format(sim_time, bow_shock_standoff))
