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
Reads all clustering.csv files in the folder and plots the average mixing 
index over time for all simulations in one graph (average is over multiple
realisations of the same simulation). 
"""

if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output_23_01/dCIL_persistence_initial_blob/A_0.5_Pe_40"
    

# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 


# Plotting style
lineStyles = ["solid", "dotted", "dashed", (0,(5,10)),"dashdot", (0,(1,1,1,1,5))]
lineColors = ["black", "blue", "orange", "green", "red", "purple", "brown", "gray", "olive", "cyan", "pink"]
styleIdx = 0
colorIdx = 0
A_dic = {}
Pe_dic = {}

# Allowed parameters: This way only certain simulations enter the graph
restrict_parameters = False # If False A_allowed is ignored
A_allowed = [0.3, 0.5, 0.7]
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
values = []


for fileName in sorted(os.listdir(FOLDER)): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+fileName[:-11]


    if paramFile is not None and paramFile[-10:]!="velocities":
        INPUT_FILE = paramFile+"_mixing.csv"
        # Find parameters
        exp = experiment(paramFile)
        params = exp.read_parameter_file(paramFile)
        A = params["areaFraction"]
        Pe = params["Pe"]

        if (A in A_allowed and Pe < 200) or restrict_parameters==False:

            # Update linestyle and line color
            if A not in A_dic:
                A_dic[A]=lineStyles[styleIdx]
                styleIdx += 1

            if Pe not in Pe_dic:
                Pe_dic[Pe]= lineColors[colorIdx]
                colorIdx += 1

            # Read mixing index from file
            mixingIdx, mixingIdxGreen, mixingIdxRed, measurementTimes = read_mixing_index_file(INPUT_FILE)

            # Plot mixing index over time
            ax_all.plot(measurementTimes, mixingIdx, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
            ax_green.plot(measurementTimes, mixingIdxGreen, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
            ax_red.plot(measurementTimes, mixingIdxRed, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])

            # Plot mixing index overview in last timestep
            tIdx = -1
            areaFractions.append(A)
            Pes.append(Pe)
            values.append(mixingIdx[-1])
        
ax_all.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) # Make room for legend on the left
ax_all.set_xlabel("Time")
ax_all.set_ylabel("Demixing index")
ax_all.legend()
fig_time_dependence.savefig(FOLDER+NAME+".pdf", format='pdf')

ax_green.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) # Make room for legend on the left
ax_green.set_xlabel("Time")
ax_green.set_ylabel("Demixing index")
ax_green.legend()
fig_time_dependence_green.savefig(FOLDER+NAME+"_green.pdf", format='pdf')

ax_red.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) # Make room for legend on the left
ax_red.set_xlabel("Time")
ax_red.set_ylabel("Demixing index")
ax_red.legend()
fig_time_dependence_red.savefig(FOLDER+NAME+"_red.pdf", format='pdf')

im = ax2.scatter(areaFractions, Pes, s=200, c=values, cmap='Blues')
ax2.set_xlabel("Area Fraction")
ax2.set_ylabel("Peclet")
fig_phasediagramm.colorbar(im)
# fig_phasediagramm.savefig(FOLDER+"mixing_phasediagram.pdf")

