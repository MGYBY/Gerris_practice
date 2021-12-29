#!/usr/bin/env python

import numpy as np

obs_loc = 120.0
x_beginning = 0.0
lx = 162.0
num_box = 54
l_box = lx/num_box
flag = 1
flag_obs = True
with open('meshFile','w') as file:
    file.write('\n')
    # first box
    bc_string = """GfsBox {
pid = 0
left = Boundary {
    BcDirichlet P depth_bc(t)
    BcDirichlet U velocity_bc(t)
    BcDirichlet V 0.0
}
top = Boundary
bottom = Boundary
}"""
    file.write(bc_string)
    file.write('\n')
    file.write("# {0}th box above".format(1))
    file.write('\n')

    # middle boxes, processor 1
    bc_string_middle1 = """GfsBox {
pid = 0
top = Boundary
bottom = Boundary
}"""

# middle boxes, processor 2
    bc_string_middle2 = """GfsBox {
pid = 1
top = Boundary
bottom = Boundary
}"""

# middle boxes, processor 3
    bc_string_middle3 = """GfsBox {
pid = 2
top = Boundary
bottom = Boundary
}"""

    # for i in range(2,27):
    #     file.write(bc_string_middle1)
    #     file.write('\n')
    #     file.write("# {0}th box above".format(i))
    #     file.write('\n')


    # for i in range(27,num_box):
    #     file.write(bc_string_middle2)
    #     file.write('\n')
    #     file.write("# {0}th box above".format(i))
    #     file.write('\n')
# the box containing obstacle, processor 4
    bc_string_middle4 = """GfsBox {
pid = 3
top = Boundary
bottom = Boundary
}"""

    # staggered pid
    for i in range(2,num_box):
        x_beginning = x_beginning+l_box
        if (x_beginning<40.0):
            file.write(bc_string_middle1)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        elif (x_beginning>=40.0 and x_beginning<80.0):
            file.write(bc_string_middle2)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        elif (x_beginning>=80.0 and x_beginning<140.0):
            file.write(bc_string_middle3)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        else:
            file.write(bc_string_middle4)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')


        # the last box still needed to be taken care of due to different BC
        # the last box
    bc_string_last = """GfsBox {
pid = 3
right = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
top = Boundary
bottom = Boundary
}"""
    file.write(bc_string_last)
    file.write('\n')
    file.write("# last box above")
    file.write('\n')

    file.write('\n')

    # box connectivity
    for i in range(1, (num_box)):
        i1 = int(i)
        i2 = int(i+1)
        file.write("{0} {1} right".format(i1, i2))
        file.write('\n')
