#!/usr/bin/env python

import numpy as np
# import time
import sys
# import re
import fnmatch
import os

normal_depth = 0.00798
obs_loc = 160.0
sim_time = sys.argv[1]
# print('Simulation time is:{0}'.format(sim_time))
# a two-column file
# a = np.loadtxt("centerlineSlice", usecols = (1,4), skiprows=1, dtype=np.float32)
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, 'centerlineSlice*'):
        a = np.loadtxt(file, usecols = (1,4), skiprows=1, dtype=np.float32)

val = a[:,1]
loc = a[:,0]
max_depth = np.max(val)
bow_shock_depth = (max_depth+normal_depth)/2.0
bow_shock_diff = np.abs(val-bow_shock_depth)
min_diff_loc = np.argmin(bow_shock_diff)
bow_shock_standoff = obs_loc-loc[min_diff_loc]
# print('Location of bow shock is:{0}'.format(bow_shock_loc))
with open('bowShockWidth','a') as file:
    file.write("{0} {1}\n".format(sim_time, bow_shock_standoff))