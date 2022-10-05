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

"""
Finds the clustering for all realisations of a given simulation and writes 
it into a file named "*_clustering.csv"
"""

# Parameter
NEIGHBOUR_DISTANCE=1.05

if(len(sys.argv)==2): #Parameter file given externally:
    paramFile = sys.argv[1]

else: #Give parameter file manually, loop over all experiments in the folder
    paramFile = "/home/marius/PhD/CellMotility/agent_simulation/new_output/test/test_file_iteration/new_test_parameters"

OUTPUT_FILE = paramFile+"_clustering.csv"
try:
    os.remove(OUTPUT_FILE)
except:
    print("")

# Calculate mixing index for every time step and realisations and save it in OUTPUT_FILE
write_data_for_all_realisations("clustering", paramFile, OUTPUT_FILE, NEIGHBOUR_DISTANCE)
