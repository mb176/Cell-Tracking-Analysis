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
from analysis_library import *
plt.style.use('seaborn-whitegrid')
plt.close("all")

"""
Finds the mixing index for all realisations of a given simulation and writes 
it into a file named "*_mixing.csv".
"""

# Parameter
NEIGHBOUR_DISTANCE=1.1

if(len(sys.argv)==2): #Parameter file given externally:
    paramFile = sys.argv[1]

else: #Give parameter file manually, loop over all experiments in the folder
    paramFile = "/home/marius/PhD/CellMotility/agent_simulation/output_23_06/CIL_based_demixing/no_cooldown_no_persistence/vary_D_A_k_5000/A_0.8_D_10"

OUTPUT_FILE = paramFile+f"_mixing_d_{NEIGHBOUR_DISTANCE}.csv"
try:
    os.remove(OUTPUT_FILE)
except:
    print("")
# Calculate mixing index for every time step and realisations and save it in OUTPUT_FILE
write_data_for_all_realisations("mixingIdx", paramFile, OUTPUT_FILE, NEIGHBOUR_DISTANCE)

print("Mixing indices written to file for all realisations of this simulation\n")
