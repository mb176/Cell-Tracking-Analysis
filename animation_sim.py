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
Creates an movie showing particles at each time step for a 
given simulation.
"""


if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output_23_01/diffCIL_persistence/movies/A_0.7_Pe_200"]
    # PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity/LowDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity/HighDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl/HighDensityControl"]


for PATH in PATHS:
    if(exists(PATH+f"_velocities_tracks_1.csv")):
        animate_tracks(PATH, velocityArrows=True)
    else:
        animate_tracks(PATH, velocityArrows=False)
 
    

