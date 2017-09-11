#!/bin/env python

import os
import numpy as np
import time



#######  SU(3) ######################################
N = 3

t = time.process_time()

# Strong coupling #####################
#for b in np.linspace(0.2,2,10):
#    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
#        " -beta=" + str(b) + \
#        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
#    print(cmd)
#    os.system(cmd)

# Weak coupling #######################
# Change range for each N
#for b in np.arange(4,20,0.5):
#    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
#        " -beta=" + str(b) + \
#        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
#    print(cmd)
#    os.system(cmd)



# beta loop ##########################
#for b in np.arange(5.0,6.5,0.1):
#    cmd = "mpirun -np 32 -machinefile procs32 su3dev2 -beta " + str(b) + \
#        " -L 8x8x8x4 -warms 5000 -trajecs 10000 -meas 5 >> out.su3dev2_l84bloop_mpi10k"
#    print(cmd)
#    os.system(cmd)


# beta loop ##########################
#for b in np.arange(3.5, 5.5, 0.25):
#    cmd = "mpirun -np 32 -machinefile procs32 su3trPaniso -beta " + str(b) + \
#        " -gamma 3.0 -L 8x8x8x8x8 -warms 2000 -trajecs 1000 -meas 5 >> out.5dL8g3bloop"
#    print(cmd)
#    os.system(cmd)


# beta loop ##########################
for b in np.arange(3.0,7.0,0.25):
    for h in np.arange(0.0, 0.1,0.02):

        cmd = "mpirun -np 32 -machinefile procs32 su3trPaniso -beta " + str(b) + \
              " -H " + str(h) + " -gamma 2.0 -L 8x8x8x8x4 -hot -warms 2500" + \
              " -trajecs 1000 -meas 5 >> out.5dL84g20bhloop"

        print(cmd)
        os.system(cmd)



#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))
