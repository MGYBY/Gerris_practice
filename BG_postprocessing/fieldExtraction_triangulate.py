import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import rc
import matplotlib
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r"\boldmath"
#plt.rcParams['text.latex.unicode'] = True


def gettingfield(fieldstr):
    exe = "gfs2oogl2D -p cartgrid.dat -c %s < %s" % (fieldstr, place)
    temp1 = sp.check_output(exe, shell=True)
    temp2 = temp1.decode("utf-8")
    temp3 = temp2.split("\n")
    fieldtemp = []
    Xtemp = []
    Ytemp = []
    for n1 in range(0, len(temp3) - 1):
        temp4 = temp3[n1].split(" ")
        Xtemp.append(float(temp4[0]))
        Ytemp.append(float(temp4[1]))
        fieldtemp.append(float(temp4[3]))
    X = np.asarray(Xtemp)
    Y = np.asarray(Ytemp)
    field = np.asarray(fieldtemp)
    X.resize((ny,nx))
    Y.resize((ny,nx))
    field.resize((ny, nx))
    return X, Y, field

def gettingfield_mod(fieldstr):
    exe = "gfs2oogl2D -p cartgrid.dat -c %s < %s" % (fieldstr, place)
    temp1 = sp.check_output(exe, shell=True)
    temp2 = temp1.decode("utf-8")
    temp3 = temp2.split("\n")
    fieldtemp = []
    Xtemp = []
    Ytemp = []
    for n1 in range(0, len(temp3) - 1):
        temp4 = temp3[n1].split(" ")
        Xtemp.append(float(temp4[0]))
        Ytemp.append(float(temp4[1]))
        fieldtemp.append(float(temp4[3]))
    X = np.asarray(Xtemp)
    Y = np.asarray(Ytemp)
    field = np.asarray(fieldtemp)
    return X, Y, field


# ----------------------------------------------------------------------------------------------------------------------


nGFS = 1
tmp = np.ones(nGFS)

folder = 'PP'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

d = 1.0
Ldomain = 5.0*d
grid_size = Ldomain/512
xmin = (-1.0)*Ldomain/2.0 - grid_size/1.0
xmax = Ldomain/2.0 + grid_size/1.0
ymin = (-1.0)*Ldomain/2.0- grid_size/1.0
ymax = Ldomain/2.0+ grid_size/1.0
nx = 512
ny = 512
x = np.linspace(xmin, xmax, num=nx)
y = np.linspace(ymin, ymax, num=ny)
z = 0
gridfile = 'cartgrid.dat'
print('saving the grid')
f = open('cartgrid.dat', 'w+')
for i in range(ny):
 for j in range(nx):
     f.write("%f %f %f\n" % (x[j], y[i], 0))
f.close()
for ti in range(nGFS):
    t = 0.2
    place = "./sim-%3.1f.gfs" % t
    if not os.path.exists(place):
        print("File not found!")
    else:
        X, Y, f = gettingfield('P')
        X1, Y1, f1 = gettingfield_mod('P')
        Xp, Yp, fp = gettingfield('P')
        X.transpose()
        Y.transpose()
        f.transpose()
        Xp.transpose()
        Yp.transpose()
        fp.transpose()
        print("Size %d,%d" % (np.size(f), np.size(fp)))
        name = "%s/%4.4d.png" %(folder, ti)
        plt.figure(figsize=(15.0, 15.0))
        rc('axes', linewidth=2)
        #plt.contour(Y/d, X/d, f, levels=3.0, colors='b')
        #plt.contour(Y/d, X/d, f, levels=[1.0, 2.0, 3.0, 4.0, 5.0], colors='b')
        plt.contourf(Xp/d, Yp/d, fp)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlabel(r'$X$', fontsize=30)
        plt.ylabel(r'$Y$', fontsize=30)
        plt.axis('square')
        plt.ylim(xmax/d, xmin/d)
        plt.xlim(ymin/d, ymax/d)
        #plt.show()
        plt.savefig(name)
        plt.close()
        print(("Done %d of %d" % (ti+1, nGFS)))

        name = "%s/tri-%4.4d.png" %(folder, ti)
        fig, ax = plt.subplots()
        #ax.set_yticks(np.arange(0,(3+1),1))
        ax.set_aspect(1)
        triang = tri.Triangulation(X1, Y1)
        tcf = ax.tricontourf(triang, f1, levels=20, vmin=0.0, cmap='bwr')
        plt.savefig(name)
        plt.close()
