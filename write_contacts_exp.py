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


def theta(phi1, phi2):
    if not math.isnan(phi1) and not math.isnan(phi2):
        theta = phi1-phi2
        if theta<0: theta+= 2*np.pi
    else: theta = None
    return theta


name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
                '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
                ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
                '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks'],
                ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
                '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']
                ]
pair = name_pairs[1]
sourceFolder = '/home/marius/PhD/CellMotility/tracking_23_01'#'/home/marius/PhD/CellMotility/tracking_ignacio_2022/'
outputFolder = '/home/marius/PhD/CellMotility/Plots/Plots_2023_01/'
subfolder = pair[0][:pair[0].rindex("/")]+"/d_9_partial_tracks/"

min_length = 0 # Min number of timesteps a track needs to have to be included
contact_radius = 9 # Matches 30 micrometer distance if 1 pixel = 2 microm
deltaT = 5 # Number of time steps that we consider to determine "instantaneous" velocity
only_full_duration_tracks= True # Exclude reference tracks that are not there for the entire simulation

# Load green tracks
green_name = pair[0][:pair[0].find('f')]
green = experiment(green_name)
green.read_xml(sourceFolder+pair[0]+'.xml', min_length)


#Load red tracks
red_name = pair[1][:pair[1].find('f')]
red = experiment(red_name)
red.read_xml(sourceFolder+pair[1]+'.xml', min_length)


experiments = [green, red] # Order matters! 
nGreenTracks = len(green.tracks)
nRedTracks = len(red.tracks)

green_t_max = green.find_t_max()
red_t_max = red.find_t_max()
assert(green_t_max==red_t_max)
t_max = green_t_max
nTracks = 0 # Keeps track of how many 

contact_ended = False
with open(outputFolder+subfolder+f"/contacts_deltaT_{deltaT}_d_{contact_radius}_log.csv", 'w') as log:
    with open(outputFolder+subfolder+f"/contacts_deltaT_{deltaT}_d_{contact_radius}.csv", 'w') as f:
        f.write(f"Duration, theta1=phiR-phi2 initial,  theta1 end, phiR_init-phiR_end, displacement_end, type, max number of contacts, contact during displacement; \n")
        f.write(f"contact_radius, {contact_radius}, deltaT, {deltaT} \n")
        # Loop over all reference tracks
        for ref_expIdx, ref_exp in enumerate(experiments):
            for ref_trackIdx, ref_track in enumerate(ref_exp.tracks):
                if len(ref_track)-1 == t_max or only_full_duration_tracks is False:
                    nTracks+=1
                    print(f"ReferenceTrack: {ref_expIdx*nGreenTracks+ref_trackIdx}\n")
                    theta1_init = {} # phiR - phi1
                    phiR_init = {} # phi1 - phi2
                    maxContacts = {} # What was the maximum number of contacts that was had during this contact?
                    contactDuration = np.zeros(nGreenTracks+nRedTracks) # Keeps track if and for how long a contact to the reference particle exists
                    last_contact_time = -np.inf # Keeps track of wether you meet another cell in the deltaT period after contacts
                    for ref_timeIdx, time in enumerate(np.array(ref_track,dtype=int)[:,0]):
                        ref_position = np.array(ref_track[ref_timeIdx][1:])
                        # Loop over all possible contacts for the reference cell
                        for expIdx, exp in enumerate(experiments):
                            for trackIdx, track in enumerate(exp.tracks):
                                # if len(track)-1 == t_max or only_full_duration_tracks is False:
                                if (ref_expIdx, ref_trackIdx)!=(expIdx, trackIdx):
                                    key = f"{expIdx},{trackIdx}" # To denote this cell in dictionaries
                                    cellIdx = expIdx*nGreenTracks+trackIdx
                                    timeIdx, = np.where(np.array(track)[:,0]==time)
                                    if len(timeIdx)>0: # Does this track exist at this time?
                                        timeIdx = timeIdx[0]
                                        position = np.array(track[timeIdx][1:])
                                        distance, phiR = exp.distance(ref_position-position, angle=True)
                                    else:
                                        distance = float('Inf')
                                        phiR = float('NaN')

                                    # Cells in contact?
                                    if distance <= contact_radius: 
                                        # New contact?
                                        if contactDuration[cellIdx] == 0: 
                                            print("New contact:", key, cellIdx, f"{ref_timeIdx}/{len(ref_track)-1}")
                                            # tmp, phi2 = exp.instanteneous_displacement(trackIdx, timeIdx, deltaT, angle=True)
                                            tmp, phi1 = ref_exp.instanteneous_displacement(ref_trackIdx, ref_timeIdx, deltaT, angle=True)
                                            theta1_init[key] = theta(phiR, phi1)
                                            phiR_init[key] = phiR
                                            maxContacts[key] = len(np.where(contactDuration!=0)[0])
                                            

                                        contactDuration[cellIdx]+=1
                                        maxContacts[key] = max(maxContacts[key],len(np.where(contactDuration!=0)[0]))

                                        if(ref_timeIdx==len(ref_track)-1): # Is the reference track ending?
                                            contact_ended = True

                                    # Did the cells use to be in contact and broke it?
                                    elif contactDuration[cellIdx] > 0: 
                                        contact_ended = True

                                    # Is the previous contacts displacement measurement disrupted by another contact?
                                    if last_contact_time >= time - deltaT and len(np.where(contactDuration!=0)[0])>0: 
                                        f.write(f", 1 \n ") # Add last entry to previous contact
                                        last_contact_time = - np.inf

                                    if contact_ended:
                                        print("Old Contact:", key, cellIdx, f"{ref_timeIdx}/{len(ref_track)-1}")
                                        maxContacts[key] = max(maxContacts[key],len(np.where(contactDuration!=0)[0]))
                                        displacement_end, phi1 = ref_exp.instanteneous_displacement(ref_trackIdx, ref_timeIdx+deltaT, deltaT, angle=True)    
                                        # tmp, phi2 = exp.instanteneous_displacement(trackIdx, timeIdx, deltaT, angle=True)
                                        theta1_end = theta(phiR, phi1)
                                        deltaPhiR = theta(phiR_init[key], phiR)
                                        # theta2 = theta(phi1, phi2)

                                        # Contact type: 0 - red&red 1 - green&green 2 - green reference with red, 3 - red reference with green
                                        if(ref_expIdx==1 and expIdx == 1): 
                                            contactType = 0
                                        elif(ref_expIdx==0 and expIdx == 0):
                                            contactType = 1
                                        elif ref_expIdx==0 and expIdx == 1:
                                            contactType = 2
                                        elif ref_expIdx==1 and expIdx == 0:
                                            contactType = 3
                                        else: 
                                            assert(False)
                                        
                                        # Write down contact information
                                        duration = contactDuration[cellIdx]
                                        print(f"{duration}, {theta1_init.get(key)}, {theta1_end}, {deltaPhiR}, {displacement_end}, {contactType}, {maxContacts.get(key)} \n")
                                        f.write(f"{duration}, {theta1_init.pop(key, None)}, {theta1_end}, {deltaPhiR}, {displacement_end}, {contactType}, {maxContacts.pop(key, None)} ")
                                        # log.write(f"Contact Start: ({ref_track[int(ref_timeIdx-duration)]}, {track[int(timeIdx-duration)]})\n")
                                        # log.write(f"Contact end: ({ref_track[int(ref_timeIdx)]}, {track[int(timeIdx)]})\n")
                                        log.write(f"Phi1:{phi1},  displacement: {displacement_end}\n")

                                        contactDuration[cellIdx] = 0
                                        contact_ended = False
                                        last_contact_time = time
                                            
                        if last_contact_time==time-deltaT or (ref_timeIdx==len(ref_track)-1 and last_contact_time>0): # Contact displacement measurement is over
                            f.write(f", 0 \n") # Add last entry to previous contact
                            last_contact_time = -np.inf



    log.write(f"Number of considered tracks: {nTracks}/{nGreenTracks+nRedTracks}")