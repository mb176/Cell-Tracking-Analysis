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
index at the final time step as a function of the area fraction.
"""

if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output_new_parameters/turnAround_persistence_t_100/areaFraction_0.7_Pe_40_mixing.csv"

# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 


# Pick parameters
Pes = [40, 80, 120, 160]
areaFractions = [0.7, 0.8]

NAME = "crossover_mixing_Pe"


# Initialise graphics
# Mixing index over time for all, green and red particles
fig, ax = plt.subplots(1,1) 

for Pe in Pes:
    mixing_values = []
    for A in areaFractions:
        INPUT_FILE = FOLDER+f"areaFraction_{A}_Pe_{Pe}_mixing.csv"

        mixingIdx, mixingIdxGreen, mixingIdxRed, measurementTimes = read_mixing_index_file(INPUT_FILE)
        mixing_values.append(mixingIdx[-1])
    
    ax.plot(areaFractions, mixing_values, label=f"Pe={Pe}")

ax.set_xlabel("Area fraction")
ax.set_ylabel("Demixing index")
ax.legend()

fig.savefig(FOLDER+NAME+".png", format='png',dpi=300)


