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
import statistics

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

min_duration = 1 # Contacts shorter or equal than this are excluded from the analysis
max_duration = 100
max_contacts = 1 # If particles had more contacts than this during the contact in question, it is excluded
CIL_angle = np.pi/6 # If the outgoing angle is smaller than this it counts as CIL
no_contacts_during_displacement = True

# Data structures
duration = np.zeros(4)
duration2 = np.zeros(4)
durationCounter = np.zeros(4)
theta1_end_angles = [[],[],[],[]]
theta1_end_counter = np.zeros(4)
CIL_counter = np.zeros(4)
displacement_end = [[],[],[],[]]
contacts_during_displacements = 0



def angle_average(angles):
    direction = np.array([np.sum(np.cos(angles)), np.sum(np.sin(angles))])
    if direction[0]>0:
        angle_avr = np.arctan(direction[1]/direction[0])
    else:
        angle_avr= np.arctan(direction[1]/direction[0]) + np.pi
    return angle_avr

with open(outputFolder+subfolder+"/contacts_deltaT_5_selected_new.csv", 'r') as f:
    csvFile = csv.reader(f, delimiter=",")
    for lineIdx, line in enumerate(csvFile):
        if lineIdx == 1:
            contact_radius = float(line[1])
            deltaT = float(line[3])
        elif lineIdx > 1 and float(line[0])>min_duration and float(line[0])<max_duration and float(line[6])<=max_contacts:
                # Contact type: 0 - red&red 1 - green&green 2 - green reference with red, 3 - red reference with green
                # "Duration, theta1=phiR-phi2 initial,  theta1 end, phiR_init-phiR_end, displacement_end, type, max number of contacts, contact during displacement; 
                #Average durations
                duration[int(float(line[5]))] += float(line[0])
                duration2[int(float(line[5]))] += float(line[0])**2
                durationCounter[int(float(line[5]))] += 1
                contacts_during_displacements = float(line[7])
                
                if line[2]!=" None":
                    if(contacts_during_displacements>0):
                        print(f"Point excluded: Type={float(line[5])}, angle={180/np.pi*float(line[2])}")
                    if(contacts_during_displacements==0 or not no_contacts_during_displacement):
                        theta1_end_angles[int(float(line[5]))].append(float(line[2]))
                        theta1_end_counter[int(float(line[5]))] += 1
                        angle = float(line[2])
                        distance = float(line[4])
                        displacement_end[int(float(line[5]))].append([np.cos(angle)*distance, np.sin(angle)*distance])
                        assert(angle>=0)
                        if(angle<CIL_angle or angle > 2*np.pi-np.pi/6):
                            CIL_counter[int(float(line[5]))] +=1
                # Average directions after contact

duration = duration/durationCounter
duration_err = np.sqrt((duration2/durationCounter - duration**2)/durationCounter)
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

# # Plot velocity angle after collision
# theta1_end_avr = []
# for typeIdx in range(4):
#     theta1_end_avr.append(360/2/np.pi*angle_average(np.array(theta1_end_angles[typeIdx])))
# fig, ax = plt.subplots()
# ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
# ax.bar([0,1,2,3], theta1_end_avr, tick_label=[f"red hom ({theta1_end_counter[0]})", f"green hom ({theta1_end_counter[1]})", f"green het ({theta1_end_counter[2]})", f"red het ({theta1_end_counter[3]})"])
# ax.set_ylabel("Angle")
# plt.xticks(rotation = 45)
# plt.tight_layout() # so labels don't get cut off
# fig.savefig(outputFolder+subfolder+"/theta1_after_contact.png",dpi=500)


# Plot CIL_counter
CIL_percentage = CIL_counter/theta1_end_counter
fig, ax = plt.subplots()
ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
ax.bar([0,1,2,3], CIL_percentage, tick_label=[f"red hom ({theta1_end_counter[0]})", f"green hom ({theta1_end_counter[1]})", f"green het ({theta1_end_counter[2]})", f"red het ({theta1_end_counter[3]})"])
ax.set_ylabel("CIL percentage")
plt.xticks(rotation = 45)
plt.tight_layout() # so labels don't get cut off
fig.savefig(outputFolder+subfolder+"/CIL_percentage.png",dpi=500)


# # Scatter plot displacement after collision
# fig, ax = plt.subplots()
# ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
# labels = ["red hom", "green hom", "green het", "red het"]
# for typeIdx in range(4):
#     points = np.array(displacement_end[typeIdx]).transpose()
#     X = points[0]
#     Y = points[1]
#     ax.scatter(X,Y, label=labels[typeIdx], s=3)
# ax.plot([-15,0],[0,0],color="black", linestyle="--",label="Last contact")
# ax.set_xlabel("X [pixel]")
# ax.set_ylabel("Y [pixel]")
# ax.legend()
# fig.savefig(outputFolder+subfolder+"/displacments_scatter.png",dpi=500)

# Plot average cumulated displacement after collision
fig, ax = plt.subplots()
ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
labels = ["red hom", "green hom", "green het", "red het"]
for typeIdx in range(4):
    points = np.array(displacement_end[typeIdx]).transpose()
    X = points[0]
    Y = points[1]
    ax.plot([0,np.sum(X)/len(X)],[0,np.sum(Y)/len(Y)],label=labels[typeIdx])
ax.plot([-15/len(X)*2,0],[0,0],color="black", linestyle="--",label="Last contact")
ax.set_xlabel("X [pixel]")
ax.set_ylabel("Y [pixel]")
plt.gca().set_aspect('equal')
ax.legend()
fig.savefig(outputFolder+subfolder+"/sum_of_displacments.png",dpi=500)

# Plot average absolute displacement
fig, ax = plt.subplots()
ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
labels = ["red hom", "green hom", "green het", "red het"]
avr_displacement = np.zeros(4)
for typeIdx in range(4):
    points = np.array(displacement_end[typeIdx])
    total_displacement = 0
    for point in points:
        total_displacement += np.linalg.norm(point)
    avr_displacement[typeIdx] = total_displacement/theta1_end_counter[typeIdx]
ax.bar([0,1,2,3], avr_displacement, tick_label=[f"red hom ({durationCounter[0]})", f"green hom ({durationCounter[1]})", f"green het ({durationCounter[2]})", f"red het ({durationCounter[3]})"])
ax.set_ylabel("Average displacment [pixel]")
plt.xticks(rotation = 45)
plt.tight_layout() # so labels don't get cut off
ax.legend()
fig.savefig(outputFolder+subfolder+f"/abs_displacement_after_contact_deltaT_{deltaT}.png",dpi=500)

# plot angle histogram of distance of theta1 from 0
fig, ax = plt.subplots(2,2)

theta1_end_abs = []
for typeIdx in range(4):
    theta1 = np.array(theta1_end_angles[typeIdx])
    theta1_prime = abs(np.array(theta1_end_angles[typeIdx])-2*np.pi)
    theta1 = np.where(theta1<np.pi,theta1, theta1_prime)
    theta1_end_abs.append(theta1)
indices = [(0,0), (0,1), (1,1), (1,0)]
labels = ["red hom", "green hom", "green het", "red het"]
for typeIdx in range(4):
    ax[indices[typeIdx]].hist(180/np.pi*np.array(theta1_end_abs[typeIdx]),label=labels[typeIdx],bins=int(np.pi/CIL_angle*2))
    ax[indices[typeIdx]].set_ylabel("Frequency")
    ax[indices[typeIdx]].legend()
# fig.savefig(outputFolder+subfolder+f"/theta1_histograms.png",dpi=500)

# Plot average angle deviation from zero
theta1_end_deviation = np.zeros(4)
for typeIdx in range(4):
    # for angle in theta1_end_angles[typeIdx]:
    #     theta1_end_deviation[typeIdx] += min(abs(angle), abs(angle-2*np.pi))
    theta1_end_deviation[typeIdx] = statistics.median(theta1_end_abs[typeIdx]) #np.sum(theta1_end_abs[typeIdx]) / theta1_end_counter[typeIdx]
fig, ax = plt.subplots()
ax.set_title(f"Min contact frames = {min_duration}, Max contacts = {max_contacts}, (R = {contact_radius}, delta T = {deltaT})")
ax.bar([0,1,2,3], 360/2/np.pi*theta1_end_deviation, tick_label=[f"red hom ({theta1_end_counter[0]})", f"green hom ({theta1_end_counter[1]})", f"green het ({theta1_end_counter[2]})", f"red het ({theta1_end_counter[3]})"])
ax.set_ylabel("Median deviation of theta1 from 0 [degrees]")
plt.xticks(rotation = 45)
plt.tight_layout() # so labels don't get cut off
# fig.savefig(outputFolder+subfolder+"/theta1_avr_deviation.png",dpi=500)



# Histograms averaged over cell types

fig, ax = plt.subplots(1,2)
theta1_hom = np.concatenate((theta1_end_abs[0], theta1_end_abs[1]), axis=0)
theta1_het = np.concatenate((theta1_end_abs[2], theta1_end_abs[3]), axis=0)
labels = ["red hom", "green hom", "green het", "red het"]
ax[0].set_title(stats.ks_2samp(theta1_hom,theta1_het)) 
ax[0].hist(180/np.pi*theta1_het,label="Heterotypic",bins=int(np.pi/CIL_angle*2))
ax[1].hist(180/np.pi*theta1_hom,label="Homotypic",bins=int(np.pi/CIL_angle*2))
ax[0].legend()
ax[1].legend()

# fig.savefig(outputFolder+subfolder+f"/theta1_histograms_2.png",dpi=500)
plt.show()