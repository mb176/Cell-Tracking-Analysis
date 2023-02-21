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
Reads all clustering.csv files in the folder and plots the average maximum
cluster size over time for all simulations in one graph (average is over multiple
realisations of the same simulation). 
"""

def get_avrMaxClusterSize(measurementTimes, clusterSizes, nParticles):
    """Plots average max cluster size for all simulations in the folder.
    Finds the average maximum cluster size at each measurement time. Finds
    the individual realisations by summing up the cluster sizes until it 
    reaches nParticles.
    """
    avrMaxClusterSize = np.zeros(len(measurementTimes))
    for tIdx in range(len(measurementTimes)):
        maxClusterSize = []
        clusters = []
        n = 0 # Keeps track of the number particles in clusters
        for size in clusterSizes[tIdx]:
            n += size
            clusters.append(size)
            if n == nParticles: # Is this realisation complete?
                maxClusterSize.append(max(clusters))
                clusters = []
                n = 0
            assert(n < nParticles), "The sum of all cluster sizes in a realisation has to be equal to the total number of particles."
        avrMaxClusterSize[tIdx] = np.sum(np.array(maxClusterSize))/len(maxClusterSize)
    return avrMaxClusterSize

if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output_23_02/persistenc_only_green_dCIL/A_0.7_Pe_40"


# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 


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
NAME = "max_cluster_log"

# Allowed parameters: This way only certain simulations enter the graph
restrict_parameters = True # If False A_allowed is ignored
A_allowed = [0.3, 0.5]


for fileName in sorted(os.listdir(FOLDER)): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+fileName[:-11]


    if paramFile is not None and paramFile[-10:]!="velocities":
        INPUT_FILE = paramFile+"_clustering.csv"
        # Find parameters
        exp = experiment(paramFile)
        params = exp.read_parameter_file(paramFile)
        A = params["areaFraction"]
        Pe = params["Pe"]
        nParticles = params["nParticles"]
        nGreenParticles = params["nGreenParticles"]
        nRedParticles = params["nRedParticles"]

        if (A in A_allowed and Pe < 200) or restrict_parameters==False:
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
            # fig_histogram.savefig(paramFile+"cluster_histogram.png", format='png',dpi=300)

            # Average maximum cluster size over all realisations
            avrMaxClusterSize = get_avrMaxClusterSize(measurementTimes, clusterSize, nParticles)
            avrMaxClusterSizeGreen = get_avrMaxClusterSize(measurementTimes, clusterSizeGreen, nGreenParticles)
            avrMaxClusterSizeRed = get_avrMaxClusterSize(measurementTimes, clusterSizeRed, nRedParticles)

            # maxClusterSize = [max(clusterSize[tIdx]) for tIdx in range(len(measurementTimes))]
            # maxClusterSizeGreen = [max(clusterSizeGreen[tIdx]) for tIdx in range(len(measurementTimes))]
            # maxClusterSizeRed = [max(clusterSizeRed[tIdx]) for tIdx in range(len(measurementTimes))]


            ax_all.plot(measurementTimes, avrMaxClusterSize, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
            ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
            ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])

# Label and save the plots
ax_all.legend(); ax_all.set_xlabel("Time"); ax_all.set_ylabel("Avr max cluster size"); ax_all.set_title("All clusters")
ax_all.set_xscale("log")
ax_all.set_yscale("log")
ax_green.legend(); ax_green.set_xlabel("Time"); ax_green.set_ylabel("Max cluster size"); ax_green.set_title("Green clusters")
ax_green.set_xscale("log")
ax_green.set_yscale("log")
ax_red.legend(); ax_red.set_xlabel("Time"); ax_red.set_ylabel("Max cluster size"); ax_red.set_title("Red clusters")
ax_red.set_xscale("log")
ax_red.set_yscale("log")

# Make room for the legend on the left
ax_all.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_green.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_red.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 


fig_all.savefig(FOLDER+NAME+".png", format='png')
fig_green.savefig(FOLDER+NAME+"_green.png", format='png')
fig_red.savefig(FOLDER+NAME+"_red.png", format='png')


