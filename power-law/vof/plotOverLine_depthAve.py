outputpersec=2
try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import numpy as np
import csv
import os

normal_depth = 0.00160556
# end point of vertical line
yy = 0.0032
hori_line_num = 400

#my_foam = FindSource("./simVTK-2.3.vtk")
#SetActiveSource(my_foam)

simVTK23vtk = LegacyVTKReader(registrationName='simVTK-2.3.vtk', FileNames=['/media/boyuan/by/gerris_runs/rw_powerLaw/vof/periodic/n04/Fr0447/minimumRW/restart/simVTK-2.3.vtk'])
renderView1 = GetActiveViewOrCreate('RenderView')

#tsteps = my_foam.TimestepValues
PlotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=simVTK23vtk, Source="Line" )
PlotOverLine1.Source.Resolution = 500
DataRepresentation7 = Show()

for xcoord in np.linspace(-0.014014, 0.084087, hori_line_num):
    PlotOverLine1.Source.Point1 = [xcoord, 0, 0]
    PlotOverLine1.Source.Point2 = [xcoord, yy, 0]
    source = PlotOverLine1
    view = GetActiveView()
    Render()
    filename = "slice_3.5.csv"
    writer = CreateWriter(filename, source)
    writer.FieldAssociation = "Point Data"
    writer.UpdatePipeline()
    Render()
    del writer
    a = np.loadtxt(filename, usecols = range(0,4), skiprows=1, dtype=np.float32, delimiter=',')
    #a = a[~np.isnan(a).any(axis=1)] # Remove the rows with NaN
    a = a[~np.isnan(a).any(axis=1), :]
    vof_array = a[:, 2]
    #vof_array = a[:, 1]
    surface_ind = 0
    # locate the free-surface
    for ele in range(np.size(vof_array)):
        if ((vof_array[ele]>0.001 and vof_array[ele]<0.999) and ele>surface_ind):
        #if (( vof_array[ele]==1 and vof_array[ele+1]==0 ) and ele>surface_ind):
            surface_ind = ele
#     sum_val = np.average((a[0:surface_ind, 2]*a[0:surface_ind, 3]))
    sum_val = np.average(a[0:surface_ind, 3], weights=vof_array[0:surface_ind])
    f=open('res','a')
    np.savetxt(f, np.array([xcoord, sum_val]), newline=" ")
    #np.savetxt(f, np.array([xcoord, sum_val, surface_ind]), newline=" ")
    #print(str(surface_ind))
    f.write("\n")
    f.close()
    os.system('rm '+ filename)
