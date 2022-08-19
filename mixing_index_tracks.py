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

# Parameters 
MIN_LENGTH = 0
NEIGHBOURS_DISTANC = 12 #two times the cell readius

file = open(TARGET_FOLDER+"/mixing_index.txt", "a")

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

    times = [0, 40, 80, 120, 160, 200]

    fig, axes = plt.subplots(1,1)

    values = []



#     #Find last time step
#     assert(green.find_t_max()==red.find_t_max()), "Tracks don't match"
#     t =green.find_t_max()

    
    for t in times:
            mixingIdx, mixing = calculate_mixing_index([green,red], t, neighbourDistance = NEIGHBOURS_DISTANC)
            values.append(mixingIdx)

    axes.plot(times,values)
    axes.set_xlabel("Time")
    axes.set_ylabel("Mixing Index")
    fig.savefig(TARGET_FOLDER+pair[0]+"_mixing.png")

#     #Get mixing index with green as reference type
#     greenMixingIdx = calculate_mixing_index([green,red], t, neighbourDistance = NEIGHBOURS_DISTANC)

#     #get it the other way around
#     redMixingIdx = calculate_mixing_index([red,green], t, neighbourDistance = NEIGHBOURS_DISTANC)

#     print(greenMixingIdx)
#     print(redMixingIdx)
#     #Write index to file
#     file.write(green_name +":   %f"%greenMixingIdx +"\n" )
#     file.write(red_name +":    %f"%redMixingIdx +"\n" )

