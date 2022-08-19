import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from analysis_library import *
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
plt.style.use('seaborn-whitegrid')
plt.close("all")

SOURCE_FOLDER="/home/marius/PhD/CellMotility/tracking_ignacio_2022"
TARGET_FOLDER="/home/marius/PhD/CellMotility/Plots/Plots_Ignacio2022"

name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
        '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
        ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
        '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks'],
        ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
        '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']
        ]

# name_pairs = [['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
#                 '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']]

#parameters
MIN_LENGTH = 0
NEIGHBOURS_DISTANCE = 16


for pair in name_pairs:
    #Load green tracks
    green_name = pair[0][:pair[0].find('f')]
    green = experiment(green_name)
    green.read_xml(SOURCE_FOLDER+pair[0]+'.xml', MIN_LENGTH)
    
    #Load red tracks
    red_name = pair[1][:pair[1].find('f')]
    red = experiment(red_name)
    red.read_xml(SOURCE_FOLDER+pair[1]+'.xml', MIN_LENGTH)
    fig, axes = plt.subplots(1,1)

    times = [0, 50, 150, 200]

    for t in times:
            fig, axes = plt.subplots(1,1)

            plot_cluster_histogram(axes, [green,red], t, neighbourDistance = NEIGHBOURS_DISTANCE)
            fig.savefig(TARGET_FOLDER+pair[0]+'cluster_histogram_d_%f_t_%f.png'%(NEIGHBOURS_DISTANCE,t), format='png',dpi=300)

    



 
