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
plt.style.use('seaborn-whitegrid')
plt.close("all")


if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/new_output/test/new_test_parameters"
            ]

#Parameters



for parameterFile in PATHS:

    # #Clustering
    # OUTPUT_FILE = parameterFile+"_clusters.csv"
    # os.remove(OUTPUT_FILE)
    # NEIGHBOUR_DISTANCE=1.05
    # write_data_for_all_realisations("clustering", parameterFile, OUTPUT_FILE, NEIGHBOUR_DISTANCE)

    # clusterSizeFrequency, measurementTimes = read_clustering_file(OUTPUT_FILE)
    # print(clusterSizeFrequency)
    
    # #Mixing
    OUTPUT_FILE = parameterFile+"_mixing.csv"
    # os.remove(OUTPUT_FILE)
    # NEIGHBOUR_DISTANCE=None
    # write_data_for_all_realisations("mixingIdx", parameterFile, OUTPUT_FILE, NEIGHBOUR_DISTANCE)

    fig, ax = fig, axes = plt.subplots(1,1)
    plot_mixing_index_over_time(ax, OUTPUT_FILE, label="Test")
    ax.legend()
    fig.savefig(parameterFile+"_mixing_idx.png", format='png',dpi=300)

    
    
    # write_clustering_for_all_realisations(PATH, OUTPUT_FILE, NEIGHBOUR_DISTANCE)
 



