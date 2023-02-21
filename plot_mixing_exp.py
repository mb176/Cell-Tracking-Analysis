import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # This avoids crashes on the math cluster from using to many resources
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
from celluloid import Camera
from scipy import spatial
import networkx as nx

from os.path import exists

# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *

"""
This file exists to test and play around with the 
functions in analysis_library.py.
"""
name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
                '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
                ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
                '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks'],
                ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
                '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']
                ]

# Select experiment to consider:
pair = name_pairs[0]
sourceFolder = '/home/marius/PhD/CellMotility/tracking_23_01'#'/home/marius/PhD/CellMotility/tracking_ignacio_2022/'
outputFolder = '/home/marius/PhD/CellMotility/Plots/Plots_2023_01'
subfolder = pair[0][:pair[0].rindex("/")]

min_length = 15
neighbourDistance = 15 # None gives voronoi tesselation

#Load green tracks
green_name = pair[0][:pair[0].find('f')]
green = experiment(green_name)
green.read_xml(sourceFolder+pair[0]+'.xml', min_length)

#Load red tracks
red_name = pair[1][:pair[1].find('f')]
red = experiment(red_name)
red.read_xml(sourceFolder+pair[1]+'.xml', min_length)

experiments = [green, red]

x_length = max(experiments[0].x_max,experiments[1].x_max)
y_length = max(experiments[0].y_max,experiments[1].y_max)
max_time = max(experiments[0].find_t_max(), experiments[1].find_t_max())

time_values = np.array(np.linspace(0,max_time,10),dtype=int)

mixingIndex = np.zeros((3,len(time_values))) # Type: Total, green, red


cellNumbers = []
for timeIdx, time in enumerate(time_values):
    mixingIndex[0,timeIdx], mixingIndex[1,timeIdx], mixingIndex[2,timeIdx]= calculate_mixing_index(experiments, time, neighbourDistance)
    tmp1, tmp2, nCells = get_positions_and_particle_types(experiments, time, count_cells=True)
    cellNumbers.append(nCells)
cellNumbers = np.array(cellNumbers)

# Plot Mixing index
fig,ax = plt.subplots(1,1)
ax.plot(2*time_values,mixingIndex[0,:], label="Overall mixing index", c="blue")
ax.plot(2*time_values,mixingIndex[1,:], label="Green mixing index", c="green")
ax.plot(2*time_values,mixingIndex[2,:], label="Red mixing index", c="red")
# ax.plot(2*time_values,cellNumbers[:,1]/cellNumbers[:,2], label="Green/red ratio", c="black", linestyle="--")
plt.legend()
ax.set_xlabel("Time [min]")
ax.set_ylabel("Demixing Index")
ax.set_title(f"Neighbour distance = {neighbourDistance}")
fig.savefig(outputFolder+subfolder+f"/mixing_d_{neighbourDistance}.png",dpi=500)

                
# Plot number of particles
fig,ax = plt.subplots(1,1)
ax.plot(2*time_values,cellNumbers[:,0], label="All cells", c="blue")
ax.plot(2*time_values,cellNumbers[:,1], label="Green cells", c="green")
ax.plot(2*time_values,cellNumbers[:,2], label="Red cells", c="red")
plt.legend()
ax.set_xlabel("Time [min]")
ax.set_ylabel("Number of cells")
ax.set_title(f"Neighbour distance = {neighbourDistance}")
fig.savefig(sourceFolder+subfolder+f"/number_of_cells.png",dpi=500)
plt.show()