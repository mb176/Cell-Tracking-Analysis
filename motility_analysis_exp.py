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
import awkward as ak
from os.path import exists
import time

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
pair = name_pairs[2]
sourceFolder = '/home/marius/PhD/CellMotility/tracking_23_01'#'/home/marius/PhD/CellMotility/tracking_ignacio_2022/'
outputFolder = '/home/marius/PhD/CellMotility/Plots/Plots_2023_01'
subfolder = pair[0][:pair[0].rindex("/")]

min_length = 31
deltaTs = [5,10,15,20,25,30]
neighbourDistance = 30 # None gives voronoi tesselation

#Load green tracks
green_name = pair[0][:pair[0].find('f')]
green = experiment(green_name)
green.read_xml(sourceFolder+pair[0]+'.xml', min_length)

#Load red tracks
red_name = pair[1][:pair[1].find('f')]
red = experiment(red_name)
red.read_xml(sourceFolder+pair[1]+'.xml', min_length)

experiments = [green, red]

#Extract parameters from experiments
x_length = max(experiments[0].x_max,experiments[1].x_max)
y_length = max(experiments[0].y_max,experiments[1].y_max)
max_time = max(experiments[0].find_t_max(), experiments[1].find_t_max())
nParticles = len(experiments[0].tracks) + len(experiments[1].tracks)
nGreenParticles = len(experiments[0].tracks)
nRedParticles = len(experiments[1].tracks)


# Find local density for each particle type
localDensities = find_local_densities(experiments, neighbourDistance)
# Exclude the missing spots in the tracks (denoted by a -1)
maskedLocalDensities = np.ma.masked_values(localDensities, -1)
avrLocalDensities = np.average(localDensities,1) #Average over all time steps (takes mask into account)
medianDensity = np.median(avrLocalDensities,0)
densityCategory = np.zeros(nParticles) # 0-low green, low red, 1- high green, low red, 2- low green, high red, 3 - high green, high red
for pIdx in range(nParticles):
    if (avrLocalDensities[pIdx,0] < medianDensity[0]) and avrLocalDensities[pIdx,1] < medianDensity[1]:
        densityCategory[pIdx] += 0
    elif (avrLocalDensities[pIdx,0] >= medianDensity[0]) and avrLocalDensities[pIdx,1] < medianDensity[1]:
        densityCategory[pIdx] += 1
    elif (avrLocalDensities[pIdx,0] < medianDensity[0]) and avrLocalDensities[pIdx,1] >= medianDensity[1]:
        densityCategory[pIdx] += 2
    elif (avrLocalDensities[pIdx,0] >= medianDensity[0]) and avrLocalDensities[pIdx,1] >= medianDensity[1]:
        densityCategory[pIdx] += 3
    
    if(pIdx >= nGreenParticles): 
        densityCategory[pIdx] += 4 # +4 for red cells



TAMSDs = []

for deltaTIdx, deltaT in enumerate(deltaTs):
    tmp, greenTAMSD, greenIndices = green.TAMSD(deltaT, min_length, overlappingIntervals=True, returnIndices=True)
    tmp, redTAMSD, redIndices = red.TAMSD(deltaT, min_length, overlappingIntervals=True, returnIndices=True)
    TAMSD = greenTAMSD + redTAMSD 
    TAMSDs.append(TAMSD)
    if deltaTIdx == 0:
        redIndices = np.array(redIndices) + nGreenParticles
        greenIndices = np.array(greenIndices)
        # trackIndices allows us to match each TAMSD with its trackIdx, and thus its densityCategory
        trackIndices = np.concatenate((greenIndices, redIndices)) 


TAMSDs = np.array(TAMSDs).transpose() # Now of the form TAMSDs[TAMSDIdx, deltaTIdx]

# Average the TAMSDs based on category 
averageTAMSD = np.zeros((len(deltaTs), 8)) # Of form avrTAMSD[deltaTIdx, categoryIdx]
counter = np.zeros((8))
for TAMSDIdx, TAMSD, in enumerate(TAMSDs):
    categoryIdx = densityCategory[trackIndices[TAMSDIdx]]
    averageTAMSD[:,categoryIdx] += TAMSD
    counter[categoryIdx] +=1

for categoryIdx in range(8):
    averageTAMSD[:,categoryIdx] = averageTAMSD[:,categoryIdx]/counter[categoryIdx]

# Fit the TAMSDs for each category

print(localDensities)


# tracks = ak.Array(concatinate())

# find_local_densities(experiments, 10)
