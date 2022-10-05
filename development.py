import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
import os
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *

"""
This file exists to test and play around with the 
functions in analysis_library.py.
"""

if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/new_output/test/test_file_iteration/new_test_parameters"
            ]

#Parameters
NEIGHBOUR_DISTANCE = None


for parameterFile in PATHS:

    # #Clustering
    OUTPUT_FILE = parameterFile+"_mixing.csv"
    tIdx = 0
    try:
        os.remove(OUTPUT_FILE)
    except:
        print("")

    write_data_for_all_realisations("mixingIdx", parameterFile, OUTPUT_FILE, NEIGHBOUR_DISTANCE)

    mixingIdx, mixingIdxGreen, mixingIdxRed, measurementTimes = read_mixing_index_file(OUTPUT_FILE)

    fig, ax = fig, axes = plt.subplots(1,1)
    # Find parameters
    exp = experiment(parameterFile)
    params = exp.read_parameter_file(parameterFile)
    A = params["areaFraction"]
    Pe = params["Pe"]
    ax.plot(measurementTimes, mixingIdx, label=f"A={A}, Pe={Pe}")
    ax.legend()
    fig.savefig(parameterFile+"_mixing.png", format='png',dpi=300)

    
    
 



