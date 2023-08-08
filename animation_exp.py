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
from celluloid import Camera


from os.path import exists

# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *

"""
This file exists to test and play around with the 
functions in analysis_library.py.
"""
name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
                '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
                ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
                '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks'],
                ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
                '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']
                ]

sourceFolder = '/home/marius/PhD/CellMotility/tracking_ignacio_2023'#'/home/marius/PhD/CellMotility/tracking_ignacio_2022/'##
# outputFolder = '/home/marius/PhD/CellMotility/Plots/Plots_2023_01'

for pair in name_pairs:
    
    subfolder = pair[0][:pair[0].rindex("/")]

    min_length = 0

    #Load green tracks
    green_name = pair[0][:pair[0].find('f')]
    green = experiment(green_name)
    green.read_xml(sourceFolder+pair[0]+'.xml', min_length)

    #Load red tracks
    red_name = pair[1][:pair[1].find('f')]
    red = experiment(red_name)
    red.read_xml(sourceFolder+pair[1]+'.xml', min_length)

    experiments = [green, red]

    x_length = max(experiments[0].x_max,experiments[1].x_max)
    y_length = max(experiments[0].y_max,experiments[1].y_max)
    max_time = max(experiments[0].find_t_max(), experiments[1].find_t_max())



    #Set up figure 
    fig,ax = plt.subplots()
    ax.axis([0, x_length, 0, y_length])
    ax.set_aspect('equal') #Makes both axis have some scaling while keeping the set limits

    R = 4.5 #Diameter of 30 micrometer, ca. 9 pixels
    transparency = 0.7 #transparency of circles 

    #Record the animation
    camera = Camera(fig)
    for time in range(max_time):
        # c = getColors(colorMeasurements[frameIdx])
        for expIdx, exp in enumerate(experiments):
            for track in exp.tracks:
                tIdx, = np.where(np.array(track)[:,0]==time)
                if len(tIdx)>0: #Does the tracked particle exist at t?
                    tIdx = tIdx[0] #Only want first (and only) occurence
                    position = track[tIdx][1:]
                    print(f"{tIdx} {position[0]} {position[1]}")
                    if expIdx==0:
                        c='green'
                    elif expIdx == 1:
                        c='red'
                    # plt.gca().plot([1,2,3],[1,4,9])
                    circle = plt.Circle((position[0],
                                        position[1]),
                                        radius=R,
                                        color=c,
                                        alpha=transparency
                                        )
                    ax.add_artist(circle)
        ax.set_title(f"Time = {max_time}")    
        camera.snap()
        #Delete all artists before next frame
        # for artist in plt.gca().lines + plt.gca().collections:
        #     artist.remove()

    animation = camera.animate()
    animation.save(sourceFolder+pair[0]+".mov")