#!/usr/bin/env python

import numpy as np

num_box = 100

with open('meshFile','w') as file:
    file.write('\n')
    # first box
    bc_string = """GfsBox {
bottom = Boundary {
    BcDirichlet U 0
}
}"""
    file.write(bc_string)
    file.write('\n')
    file.write("# {0}th box above".format(1))
    file.write('\n')

    # middle boxes, processor 1
    bc_string_middle1 = """GfsBox { }"""
    # staggered pid
    for i in range(2,num_box):
        file.write(bc_string_middle1)
        file.write('\n')
        file.write("# {0}th box above".format(i))
        file.write('\n')

        # the last box still needed to be taken care of due to different BC
        # the last box
    bc_string_last = """GfsBox {
    top = Boundary
}"""
    file.write(bc_string_last)
    file.write('\n')
    file.write("# last box above")
    file.write('\n')

    # box connectivity
    for i in range(0, (num_box)):
        i1 = int(i+1)
        i2 = int(i+1)
        file.write("{0} {1} right".format(i1, i2))
        file.write('\n')

    for i in range(1, (num_box)):
        i1 = int(i)
        i2 = int(i+1)
        file.write("{0} {1} top".format(i1, i2))
        file.write('\n')
