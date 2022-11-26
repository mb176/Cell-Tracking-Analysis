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
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
plt.close("all")

"""
Reads all clustering.csv files in the folder and plots the average
cluster size over time for all simulations in one graph.
"""

def get_avrClusterSize(measurementTimes, clusterSizes):
    """Return average cluster size for each time step.
    """
    avrClusterSize = np.zeros(len(measurementTimes))
    cutoff = 0
    for tIdx in range(len(measurementTimes)):
        sizes = [x for x in clusterSizes[tIdx] if x > cutoff]
        avrClusterSize[tIdx] = np.array(sizes).sum()/len(sizes)
    return avrClusterSize

if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/new_output/turnAround_persistence_t_1000/areaFraction_0.5_Pe_120"

# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 

# Parameter
NEIGHBOUR_DISTANCE=None
NAME = "avr_cluster"

# Initialise graphic
fig_all, ax_all = plt.subplots(1,1)
fig_green, ax_green = plt.subplots(1,1)
fig_red, ax_red = plt.subplots(1,1)

# Plotting style
lineStyles = ["solid", "dotted", "dashed", (0,(5,10)),"dashdot", (0,(1,1,1,1,5))]
lineColors = ["black", "blue", "orange", "green", "red", "purple", "brown", "gray", "olive", "cyan", "pink"]
styleIdx = 0
colorIdx = 0
A_dic = {}
Pe_dic = {}

# Allowed parameters: This way only certain simulations enter the graph
restrict_parameters = False # If False A_allowed is ignored
A_allowed = [0.4, 0.5, 0.6]


for fileName in os.listdir(FOLDER): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+fileName[:-11]


    if paramFile != None:
        INPUT_FILE = paramFile+"_clustering.csv"
        # Find parameters
        exp = experiment(paramFile)
        params = exp.read_parameter_file(paramFile)
        A = params["areaFraction"]
        Pe = params["Pe"]
        nParticles = params["nParticles"]
        nGreenParticles = params["nGreenParticles"]
        nRedParticles = params["nRedParticles"]

        if (A in A_allowed) or restrict_parameters==False:
            # Update linestyle and line color
            if A not in A_dic:
                A_dic[A]=lineStyles[styleIdx]
                styleIdx += 1

            if Pe not in Pe_dic:
                Pe_dic[Pe]= lineColors[colorIdx]
                colorIdx += 1

            # Read cluster sizes from file
            clusterSize, clusterSizeGreen, clusterSizeRed, measurementTimes = read_clustering_file(INPUT_FILE)
            # # Plot  cluster size histogram for each experiment
            # fig_histogram, ax1 = plt.subplots(1,1)
            # tIdx = -1
            # plot_cluster_histogram(ax1, clusterSizeFrequencies[tIdx], label=f"A={A}, Pe={Pe}" ,color="green")        
            # ax1.legend()
            # fig_histogram.savefig(paramFile+"cluster_histogram.pdf", format='pdf')

            # Average maximum cluster size over all realisations
            avrClusterSize = get_avrClusterSize(measurementTimes, clusterSize)
            avrClusterSizeGreen = get_avrClusterSize(measurementTimes, clusterSizeGreen)
            avrClusterSizeRed = get_avrClusterSize(measurementTimes, clusterSizeRed)

            # maxClusterSize = [max(clusterSize[tIdx]) for tIdx in range(len(measurementTimes))]
            # maxClusterSizeGreen = [max(clusterSizeGreen[tIdx]) for tIdx in range(len(measurementTimes))]
            # maxClusterSizeRed = [max(clusterSizeRed[tIdx]) for tIdx in range(len(measurementTimes))]


            ax_all.plot(measurementTimes, avrClusterSize, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
            ax_green.plot(measurementTimes, avrClusterSizeGreen, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
            ax_red.plot(measurementTimes, avrClusterSizeRed, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])

# Label and save the plots
ax_all.legend(); ax_all.set_xlabel("Time"); ax_all.set_ylabel("Average cluster size"); ax_all.set_title("All clusters")
ax_all.set_xscale("log")
ax_all.set_yscale("log")
ax_green.legend(); ax_green.set_xlabel("Time"); ax_green.set_ylabel("Average cluster size"); ax_green.set_title("Green clusters")
ax_green.set_xscale("log")
ax_green.set_yscale("log")
ax_red.legend(); ax_red.set_xlabel("Time"); ax_red.set_ylabel("Average cluster size"); ax_red.set_title("Red clusters")
ax_red.set_xscale("log")
ax_red.set_yscale("log")

# Make room for the legend on the left
ax_all.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_green.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_red.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 

# add Reference lines
reference = np.array(measurementTimes)**(1/2)
ax_all.plot(measurementTimes, reference, color="black", linestyle="dashed")
ax_green.plot(measurementTimes, reference, color="black", linestyle="dashed")
ax_red.plot(measurementTimes, reference, color="black", linestyle="dashed")


fig_all.savefig(FOLDER+NAME+".pdf", format='pdf')
fig_green.savefig(FOLDER+NAME+"_green.pdf", format='pdf')
fig_red.savefig(FOLDER+NAME+"_red.pdf", format='pdf')


