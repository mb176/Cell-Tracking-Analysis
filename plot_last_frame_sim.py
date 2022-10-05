import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # This avoids crashes on the math cluster from using to many resources
import matplotlib.pyplot as plt 
import numpy as np
import csv
import os
import matplotlib.animation as animation
from celluloid import Camera
from matplotlib.collections import PatchCollection
import sys
sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
#python celluloid
#You can use np.arrays with named vectors, just like pandas!
# rain_drops = np.zeros(n_drops, dtype=[('position', float, (2,)),
#                                       ('size',     float),
#                                       ('growth',   float),
#                                       ('color',    float, (4,))])

# # Initialize the raindrops in random positions and with
# # random growth rates.
# rain_drops['position'] = np.random.uniform(0, 1, (n_drops, 2))

"""
Plots and saves a picture of the final arrangement of particles.
"""

if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/new_output/test/test_file_iteration/new_test_parameters"]
    # PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity/LowDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity/HighDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl/HighDensityControl"]

for PATH in PATHS:
    final_snapshot(PATH)
    

