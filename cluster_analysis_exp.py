import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from analysis_library import *
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
plt.style.use('seaborn-whitegrid')
plt.close("all")

"""
Calculate the cluster sizes and plots a histogram for a given experiment.
""" 

SOURCE_FOLDER="/home/marius/PhD/CellMotility/tracking_ignacio_2022"
TARGET_FOLDER="/home/marius/PhD/CellMotility/Plots/Plots_Ignacio2022"

name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
        '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
        ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
        '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks'],
        ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
        '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']
        ]

name_pairs = [['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
                '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']]

#parameters
MIN_LENGTH = 0
NEIGHBOURS_DISTANCE = 16
times = [100] # The time



for pair in name_pairs:
    # Load green tracks
    green_name = pair[0][:pair[0].find('f')]
    green = experiment(green_name)
    green.read_xml(SOURCE_FOLDER+pair[0]+'.xml', MIN_LENGTH)
    
    # Load red tracks
    red_name = pair[1][:pair[1].find('f')]
    red = experiment(red_name)
    red.read_xml(SOURCE_FOLDER+pair[1]+'.xml', MIN_LENGTH)
    fig, axes = plt.subplots(1,1)

    # Add colors to the experiments (write_cluster_to_file expects these)
    tmp1 = [" red"]*(max(times)+1)
    red.color = [tmp1] * len(red.tracks)
    tmp2 = [" green"]*(max(times)+1)
    green.color = [tmp2] * len(red.tracks)

    SUBFOLDER = pair[0][:pair[0].rindex("/")]
    OUTPUT_FILE = TARGET_FOLDER+SUBFOLDER+"/clustering.csv"
    try: # Delete old file if it already exists
                os.remove(OUTPUT_FILE)
    except:
                print("")

    

    for timeIdx, time in enumerate(times):   
        # Find clustering and write to file:
        write_clusters_to_file(OUTPUT_FILE, [green, red], time, NEIGHBOURS_DISTANCE, newRealisation=False)

    # Read clustering from file and plot histgram
    clusterSizes, clusterSizesGreen, clusterSizesRed, measurementTimes = read_clustering_file(OUTPUT_FILE)

    for timeIdx, time in enumerate(times):   
        # Plot histgrams
        fig_histogram, ax_histogram = plt.subplots(1,1)
        histogramGreen = plot_cluster_histogram(ax_histogram, clusterSizesGreen[timeIdx], label=f"Green cells" ,color="green")   
        histogramRed = plot_cluster_histogram(ax_histogram, clusterSizesRed[timeIdx], label=f"Red cells" ,color="red")

        # Test for statistically significant differences between red and green distributions
        significance = stats.ks_2samp(histogramGreen, histogramRed)  
        ax_histogram.set_title(f"Time={measurementTimes[timeIdx]}, {significance}")
        fig_histogram.savefig(TARGET_FOLDER+SUBFOLDER+'/cluster_histogram_d_%f_t_%f.png'%(NEIGHBOURS_DISTANCE,time), format='png',dpi=300)

    



 
