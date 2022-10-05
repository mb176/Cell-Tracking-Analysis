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
Reads the clustering.csv file for a given simulation and plots the cluster 
size histogram for red and green particles. Also checks for statistically
relevant differences between the distributions.
"""

if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/new_output/turnAround_persistence/areaFraction_0.5_Pe_40"
INPUT_FILE = PATH+"_clustering.csv"


# Find parameters
exp = experiment(PATH)
params = exp.read_parameter_file(PATH)
A = params["areaFraction"]
Pe = params["Pe"]
nParticles = params["nParticles"]
nGreenParticles = params["nGreenParticles"]
nRedParticles = params["nRedParticles"]
timeIdx = -1


# Read cluster sizes from file
clusterSizes, clusterSizesGreen, clusterSizesRed, measurementTimes = read_clustering_file(INPUT_FILE)

# Plot histogram
fig_histogram, ax_histogram = plt.subplots(1,1)
# plot_cluster_histogram(ax_histogram, clusterSizes[timeIdx], label=f"All cells" ,color="blue")   
histogramGreen = plot_cluster_histogram(ax_histogram, clusterSizesGreen[timeIdx], label=f"Green cells" ,color="green")   
histogramRed = plot_cluster_histogram(ax_histogram, clusterSizesRed[timeIdx], label=f"Red cells" ,color="red")  

# Test for statistically significant differences between red and green distributions
significance = stats.ks_2samp(histogramGreen, histogramRed)

ax_histogram.set_title(f"Time={measurementTimes[timeIdx]}, {significance}")
ax_histogram.legend()
fig_histogram.savefig(PATH+"_cluster_histogram.png", format='png',dpi=300)


