#!/usr/bin/env python

import numpy as np

# 6 cores for domain decomposition
# free-slip top BC
# inlet-outlet
x_beginning = 0.0
lx = 9.75
num_box = 37
l_box = lx/num_box

num_core = 6
seg_length = lx/num_core

count = 0

with open('meshFile','w') as file:
    file.write('\n')
    # first box
    bc_string = """GfsBox {
pid = 0
left = Boundary {
    BcDirichlet P hydroPressureDist(y, t)
    BcDirichlet T depth_bc(y, t)
    BcDirichlet U velocity_bc(y, t)
    BcDirichlet V 0.0
}
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""
    file.write(bc_string)
    file.write('\n')
    file.write("# {0}th box above".format(1))
    file.write('\n')
    count += 1

    # middle boxes, processor 1
    bc_string_middle1 = """GfsBox {
pid = 0
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""

# middle boxes, processor 2
    bc_string_middle2 = """GfsBox {
pid = 1
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""

# middle boxes, processor 3
    bc_string_middle3 = """GfsBox {
pid = 2
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""

# middle boxes, processor 4
    bc_string_middle4 = """GfsBox {
pid = 3
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""

# middle boxes, processor 5
    bc_string_middle5 = """GfsBox {
pid = 4
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""

# middle boxes, processor 6
    bc_string_middle6 = """GfsBox {
pid = 5
top = Boundary {
	BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
}
}"""

    # uniformly-distributed pid, reduce communication time
    for i in range(2,num_box):
        x_beginning = x_beginning+l_box
        count += 1
        if (count<=6):
            file.write(bc_string_middle1)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        elif (count>=7 and count<12):
            file.write(bc_string_middle2)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        elif (count>=12 and count<18):
            file.write(bc_string_middle3)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        elif (count>=18 and count<24):
            file.write(bc_string_middle4)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        elif (count>=24 and count<30):
            file.write(bc_string_middle5)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')
        else:
            file.write(bc_string_middle6)
            file.write('\n')
            file.write("# {0}th box above".format(i))
            file.write('\n')


        # the last box still needed to be taken care of due to different BC
        # the last box
    bc_string_last = """GfsBox {
pid = 5
right = Boundary {
    BcNeumann P 0.0
    BcNeumann T 0.0
    BcNeumann U 0.0
    BcNeumann V 0.0
}
top = Boundary {
    BcDirichlet P (-AIRRHO*GRAVITYCOEFF*CHANNELCOS*(LBOX-NORMALDEPTH))
	# BcDirichlet P 0.0
	BcDirichlet V 0.0
}
bottom = Boundary {
    BcDirichlet U 0.0
    BcDirichlet V 0.0
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
