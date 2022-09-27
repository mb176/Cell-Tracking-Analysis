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


if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/test/exampleSystem"
            ]


#Parameters
NEIGHBOUR_DISTANCE = 1.05

for PATH in PATHS:
    green, red, params = load_simulation_data(PATH)

    # Find last time step
    assert(green.find_t_max()==red.find_t_max()), "Tracks don't finish at the same time"
    t_final = green.find_t_max()

    #Find largest cluster
    maxClusterSize = max_cluster_size([green],t_final, neighbourDistance = NEIGHBOUR_DISTANCE)
    maxClusterSize = max(maxClusterSize,max_cluster_size([red],t_final, neighbourDistance = NEIGHBOUR_DISTANCE))

    fig, axes = plt.subplots(1,1)
    histogramGreen = plot_cluster_histogram(axes, [green], t_final, "Green cells" ,"g", neighbourDistance = NEIGHBOUR_DISTANCE, nBins=maxClusterSize)
    histogramRed   = plot_cluster_histogram(axes, [red], t_final, "Red cells" ,"r", neighbourDistance = NEIGHBOUR_DISTANCE, nBins=maxClusterSize)

    #Statistical relevance
    significance = stats.ks_2samp(histogramGreen, histogramRed)
    #axes.set_title(f"t={t}, N={G.number_of_nodes()}")

    axes.set_title(significance)

    fig.savefig(PATH+'_cluster_histogram_d_%f_t_%f.png'%(NEIGHBOUR_DISTANCE,t_final), format='png',dpi=300)
    
    # #Green clusters
    # fig, axes = plt.subplots(1,1)
    # plot_cluster_histogram(axes, [green], t_final, neighbourDistance = neighbourDistance)
    # fig.savefig(PATH+'_cluster_histogram_green_d_%f_t_%f.png'%(neighbourDistance,t_final), format='png',dpi=300)

    # #red clusters
    # fig, axes = plt.subplots(1,1)
    # plot_cluster_histogram(axes, [red], t_final, neighbourDistance = neighbourDistance)
    # fig.savefig(PATH+'_cluster_histogram_red_d_%f_t_%f.png'%(neighbourDistance,t_final), format='png',dpi=300)


    plt.close('all')