import subprocess
import decimal
import os
import numpy as np
 
file_delta = 0.005
file_begin = 7.995
file_end = 8.2+file_delta
#gerris2D -e "OutputSimulation { istep = 1 } simTec-${time}.plt { format = Tecplot variables = T,Tf,U,V,P,UT,VT }" ./refinedOutput/snapshot-${time}.gfs > /dev/null
for f in np.arange(file_begin, file_end, file_delta):
    f1 = float("{:.3f}".format(f))
    time_norm = decimal.Decimal(str(f1)).normalize()
    time_normStr = str(time_norm)
    print(time_normStr)
    command_str = 'gerris2D -e "OutputSimulation { istep = 1 } simTec-'+time_normStr+'.plt { format = Tecplot variables = T,Tf,U,V,P,UT,VT }" ./refinedOutput/snapshot-'+time_normStr+'.gfs > /dev/null'
    os.system(command_str)
