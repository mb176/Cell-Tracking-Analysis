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
Reads all clustering.csv and mixing.csv files in the folder and . 
"""

if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output_23_03/phasediagram/t_100/A_0.8_Pe_240"
    

# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 



# Allowed parameters: This way only certain simulations enter the graph
restrict_parameters = True # If False A_allowed is ignored
A_allowed = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
MIXING_FILE_NAME = "_mixing_d_None.csv" #"_mixing.csv""
CLUSTERING_FILE_NAME = "_clustering.csv"
NAME = "mixing_over_time"

# Initialise graphics
# Mixing index over time for all, green and red particles
fig_time_dependence, ax_all = plt.subplots(1,1) 
fig_time_dependence_green, ax_green= plt.subplots(1,1)
fig_time_dependence_red, ax_red = plt.subplots(1,1)
# Phasediagramm of the mixing index in the last timestep for all particles
fig_phasediagramm, ax2 = plt.subplots(1,1)
areaFractions = []
Pes = []
mixing = []
clusterRed = []
clusterGreen = []



for fileName in sorted(os.listdir(FOLDER)): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+fileName[:-11]


    if paramFile is not None and paramFile[-10:]!="velocities":
        # Find parameters
        exp = experiment(paramFile)
        params = exp.read_parameter_file(paramFile)
        A = params["areaFraction"]
        Pe = params["Pe"]
        nGreenParticles=params["nGreenParticles"]
        nRedParticles=params["nRedParticles"]

        if (A in A_allowed and Pe < 300) or restrict_parameters==False:

            # Read mixing index from file
            mixingIdx, mixingIdxGreen, mixingIdxRed, measurementTimes = read_mixing_index_file(paramFile+MIXING_FILE_NAME)

            # Read clustering from file
            clusterSize, clusterSizeGreen, clusterSizeRed, measurementTimes = read_clustering_file(paramFile+CLUSTERING_FILE_NAME)
            avrMaxClusterSizeGreen = get_avrMaxClusterSize(measurementTimes, clusterSizeGreen, nGreenParticles)
            avrMaxClusterSizeRed = get_avrMaxClusterSize(measurementTimes, clusterSizeRed, nRedParticles)

            # Plot mixing index overview in last timestep
            tIdx = -1
            areaFractions.append(A)
            Pes.append(Pe)
            mixing.append(mixingIdx[tIdx])
            clusterGreen.append(avrMaxClusterSizeGreen[tIdx])
            clusterRed.append(avrMaxClusterSizeRed[tIdx])

cm = plt.cm.get_cmap('RdYlBu') # get color map
maxCluster = max(max(clusterGreen), max(clusterRed))

fig, ax = plt.subplots(1,1)
ax.set_title("Demixing index")
sc = ax.scatter(Pes, areaFractions, c=mixing, s=150, cmap=cm) #, vmin=0, vmax=20, s=35, cmap=cm
ax.set_xlabel("Pe")
ax.set_ylabel("Area Fraction")
plt.colorbar(sc)
fig.savefig(FOLDER+f"mixing_diagram.png",dpi=500)


fig, ax = plt.subplots(1,1)
ax.set_title("Size of largest green cluster")
sc = ax.scatter(Pes, areaFractions, c=clusterGreen, s=150, vmax=maxCluster, cmap=cm) #, vmin=0, vmax=20, s=35, cmap=cm
ax.set_xlabel("Pe")
ax.set_ylabel("Area Fraction")
plt.colorbar(sc)
fig.savefig(FOLDER+f"green_cluster_diagram.png",dpi=500)


fig, ax = plt.subplots(1,1)
ax.set_title("Red clusters")
sc = ax.scatter(Pes, areaFractions, c=clusterRed, s=150, vmax=maxCluster, cmap=cm) #, vmin=0, vmax=20, s=35, cmap=cm
plt.colorbar(sc)

fig, ax = plt.subplots(1,1)
ax.set_title("Size difference between largest green and red cluster")
cluster_difference = np.array(clusterGreen)-np.array(clusterRed)
sc = ax.scatter(Pes, areaFractions, c=cluster_difference, s=150, cmap=cm) #, vmin=0, vmax=20, s=35, cmap=cm
ax.set_xlabel("Pe")
ax.set_ylabel("Area Fraction")
plt.colorbar(sc)
fig.savefig(FOLDER+f"cluster_difference_diagram.png",dpi=500)

plt.show()
print("f")
        

