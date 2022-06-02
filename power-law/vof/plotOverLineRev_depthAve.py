outputpersec=2
try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import numpy as np
import csv
import os

normal_depth = 0.00160556
# end point of vertical line
xx = 0.0032
line_num = 50

#my_foam = FindSource("./simVTK-2.3.vtk")
#SetActiveSource(my_foam)

simVTK23vtk = LegacyVTKReader(registrationName='simVTK-2.3.vtk', FileNames=['/media/boyuan/by/gerris_runs/rw_powerLaw/vof/periodic/n04/Fr0447/minimumRW/restart/simVTK-2.3.vtk'])
renderView1 = GetActiveViewOrCreate('RenderView')

#tsteps = my_foam.TimestepValues
PlotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=simVTK23vtk, Source="Line" )
DataRepresentation7 = Show()

for xcoord in np.linspace(-0.014014, 0.084087, line_num):
    PlotOverLine1.Source.Point1 = [0, 0, 0]
    PlotOverLine1.Source.Point2 = [0, xx, 0]
    #select fields
    #passArrays1 =PassArrays(Input=PlotOverLine1)
    ##passArrays1.PointDataArrays = ['U','alpha.water']
    #passArrays1.PointDataArrays = ['T']
    source = PlotOverLine1

    view = GetActiveView()
    #view.ViewTime = tsteps[TimeStepNum]
    Render()
    writer = CreateWriter("slice_2.3.csv", source)
    writer.FieldAssociation = "Point Data"
    writer.UpdatePipeline()
    Render()
    del writer
    a = np.loadtxt("slice_2.3.csv", usecols = range(0,4), skiprows=1, dtype=np.float32, delimiter=',')
    sum_val = np.average(a[:, 1])
    f=open('res','a')
    np.savetxt(f, np.array([xx, sum_val]), newline=" ")
    f.write("\n")
    f.close()
    os.system('rm slice_2.3.csv')


#with open("ave_val",'w') as file:
    #print("Writing filtered data to file ... ...")
    #writer = csv.writer(file, delimiter='\t')
    #np.savetxt("fn", np.transpose(arr),newline=" ")

    #writer.writerow(str(sum_val))


##for TimeStepNum in range(0,len(tsteps)):
#for TimeStepNum in range(0,403,1):
    #view = GetActiveView()
    #view.ViewTime = tsteps[TimeStepNum]
    #Render()
    #writer = CreateWriter("%dws1.csv" %(TimeStepNum), source)
    #writer.FieldAssociation = "Point Data"
    #writer.UpdatePipeline()
    #Render()
    #del writer
