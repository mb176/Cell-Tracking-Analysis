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
from scipy import spatial
import networkx as nx
import math

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
pair = name_pairs[2]
sourceFolder = '/home/marius/PhD/CellMotility/tracking_23_01'#'/home/marius/PhD/CellMotility/tracking_ignacio_2022/'
outputFolder = '/home/marius/PhD/CellMotility/Plots/Plots_2023_01'
subfolder = pair[0][:pair[0].rindex("/")]

min_duration = 2 # Contacts shorter or equal than this are excluded from the analysis
max_contacts = 1 # If particles had more contacts than this during the contact in question, it is excluded

# Data structures
duration = np.zeros(4)
duration2 = np.zeros(4)
durationCounter = np.zeros(4)
theta1_end_angles = [[],[],[],[]]
theta1_end_counter = np.zeros(4)


def angle_average(angles):
    direction = np.array([np.sum(np.cos(angles)), np.sum(np.sin(angles))])/len(angles)
    if direction[0]>0:
        angle_avr = np.arctan(direction[1]/direction[0])
    else:
        angle_avr= np.arctan(direction[1]/direction[0]) + np.pi
    return angle_avr

with open(outputFolder+subfolder+"/contacts.csv", 'r') as f:
    csvFile = csv.reader(f, delimiter=",")
    for lineIdx, line in enumerate(csvFile):
        if lineIdx == 1:
            contact_radius = float(line[1])
            deltaT = float(line[3])
        elif lineIdx > 1 and float(line[0])>min_duration and float(line[4][3])<=max_contacts:
                # Contact type: 0 - red&red 1 - green&green 2 - green reference with red, 3 - red reference with green
                #Average durations
                duration[int(float(line[4][1]))] += float(line[0])
                duration2[int(float(line[4][1]))] += float(line[0])**2
                durationCounter[int(float(line[4][1]))] += 1
                if line[2]!=" None":
                    theta1_end_angles[int(float(line[4][1]))].append(float(line[2]))
                    theta1_end_counter[int(float(line[4][1]))] += 1
                # Average directions after contact

duration = duration/durationCounter
duration_err = np.sqrt(duration2/durationCounter - duration**2)
print(duration2/durationCounter, duration_err)

# Plot duration
fig, ax = plt.subplots()
ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
ax.bar([0,1,2,3], 2*duration, tick_label=[f"red hom ({durationCounter[0]})", f"green hom ({durationCounter[1]})", f"green het ({durationCounter[2]})", f"red het ({durationCounter[3]})"])
ax.errorbar([0,1,2,3], 2*duration, yerr=2*duration_err,fmt="o", color="r")
ax.set_ylabel("Time [min]")
plt.xticks(rotation = 45)
plt.tight_layout() # so labels don't get cut off
fig.savefig(outputFolder+subfolder+"/contact_duration.png",dpi=500)
plt.show()

# Plot velocity angle after collision
theta1_end_avr = []
for typeIdx in range(4):
    theta1_end_avr.append(360/2/np.pi*angle_average(np.array(theta1_end_angles[typeIdx])))
fig, ax = plt.subplots()
ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
ax.bar([0,1,2,3], theta1_end_avr, tick_label=[f"red hom ({theta1_end_counter[0]})", f"green hom ({theta1_end_counter[1]})", f"green het ({theta1_end_counter[2]})", f"red het ({theta1_end_counter[3]})"])
ax.set_ylabel("Angle")
plt.xticks(rotation = 45)
plt.tight_layout() # so labels don't get cut off
fig.savefig(outputFolder+subfolder+"/theta1_after_contact.png",dpi=500)
plt.show()


        