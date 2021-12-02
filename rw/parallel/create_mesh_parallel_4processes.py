#!/usr/bin/env python

import numpy as np

num_box = 54
flag = 0
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
top = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
bottom = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
}"""
    file.write(bc_string)
    file.write('\n')
    file.write("# {0}th box above".format(1))
    file.write('\n')
    flag = flag+1

    # middle boxes
    bc_string_middle1 = """GfsBox {
pid = 0
top = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
bottom = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
}"""

# middle boxes
    bc_string_middle2 = """GfsBox {
pid = 1
top = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
bottom = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
}"""

# middle boxes
    bc_string_middle3 = """GfsBox {
pid = 2
top = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
bottom = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
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

    # staggered pid
    for i in range(2,num_box):
        if (flag%3==0):
            file.write(bc_string_middle2)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
            flag = flag+1
        elif(flag%3==1):
            file.write(bc_string_middle1)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
            flag = flag+1
        else:
            file.write(bc_string_middle3)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
            flag = flag+1

    # the last box
    bc_string_last = """GfsBox {
pid = 3
right = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
top = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
bottom = Boundary {
    BcNeumann P 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
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
