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
import pandas as pd
from celluloid import Camera
import scipy as sp
from scipy import optimize
from scipy import spatial
import networkx as nx
import math

from os.path import exists

# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *


def get_contact_type(ref_expIdx, expIdx):
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
    return contactType
    
######################## Parameters ############################
contact_radius = 9 # Defines what we mean by contact
cutoff_time = 15 # Cuts off the measurement here
min_length = 0 # Exclude tracks with less than this length
min_measurement_duration = 15 # Contacts that result in free paths shorter than that are excluded
only_full_duration_tracks = False 

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

# Observables, each measurement adds to the total
distance = np.zeros((4,t_max)) #[type of contact, time idx]
displacement = np.zeros((4,t_max))
displacement_var = np.zeros((4,t_max)) # Variance of displacement
MSD = np.zeros((4,t_max))
MSD_var = np.zeros((4,t_max)) #Variance of MSD
contact_cnt = np.zeros((4,t_max)) #Keeps track of total number of contacts

# Helper variables
n_measurements = np.zeros(4)




# Loop over all reference tracks
for ref_expIdx, ref_exp in enumerate(experiments):
    for ref_trackIdx, ref_track in enumerate(ref_exp.tracks):
        n_previous_contacts = -1 # To avoid measurement starting immediately
        measuring = False
        if True:
            if len(ref_track)-1 == t_max or only_full_duration_tracks is False:
                nTracks+=1
                print(f"Track number {nTracks}")
                for ref_timeIdx, time in enumerate(np.array(ref_track,dtype=int)[:,0]):
                    ref_position = np.array(ref_track[ref_timeIdx][1:])
                    #Find out if the reference cell is in contact with another cell
                    n_contacts = 0
                    for expIdx, exp in enumerate(experiments):
                        for trackIdx, track in enumerate(exp.tracks):
                            timeIdx, = np.where(np.array(track)[:,0]==time)
                            if (ref_expIdx, ref_trackIdx)!=(expIdx, trackIdx):
                                if len(timeIdx)>0: # Does this track exist at this time?
                                    timeIdx = timeIdx[0]
                                    position = np.array(track[timeIdx][1:])
                                    d = exp.distance(ref_position-position)
                                else:
                                    d = float('Inf')
                                if d <= contact_radius: 
                                    n_contacts += 1
                                    contact_type = get_contact_type(ref_expIdx, expIdx)

                    if n_contacts == 0 and n_previous_contacts==1: # Start new measurement
                        measuring = True
                        n_measurements[contact_type] += 1
                        position_at_last_contact = ref_position
                        #Reset observables
                        contact_duration = 0
                        new_displacement = 0
                        previous_position = ref_position
                        tmp_distance = np.zeros(t_max)
                        tmp_displacement = np.zeros(t_max)
                        tmp_displacement_var = np.zeros(t_max)
                        tmp_MSD = np.zeros(t_max)
                        tmp_MSD_var = np.zeros(t_max)
                        tmp_contact_cnt = np.zeros(t_max)
                        tmp_contact_cnt[contact_duration] += 1
                        old_contact_type = contact_type
                    elif n_contacts == 0 and measuring: # Continue Measurement
                        contact_duration += 1
                        new_distance = ref_exp.distance(ref_position-position_at_last_contact)
                        new_displacement += ref_exp.distance(ref_position-previous_position)
                        previous_position = ref_position

                        #Take preliminary measurements
                        tmp_distance[contact_duration] += new_distance
                        tmp_MSD[contact_duration] += new_distance**2
                        tmp_MSD_var[contact_duration] += new_distance**4
                        tmp_displacement[contact_duration] += new_displacement
                        tmp_displacement_var[contact_duration] += new_displacement**2
                        tmp_contact_cnt[contact_duration] += 1
                    elif n_contacts > 0 and measuring: # Measurement ended
                        # Save preliminary measurements if track is suitable
                        if contact_duration >= min_measurement_duration:
                            distance[old_contact_type,:] += tmp_distance
                            displacement[old_contact_type,:] += tmp_displacement
                            displacement_var[old_contact_type,:] += tmp_displacement_var
                            MSD[old_contact_type,:] += tmp_MSD
                            MSD_var[old_contact_type,:] += tmp_MSD_var
                            contact_cnt[old_contact_type,:] += tmp_contact_cnt
                    else:
                        measuring = False

                    
                    n_previous_contacts = n_contacts

# # Save results to file via pandas dataframe
# df = pd.DataFrame({ 'contact_cnt red hom': contact_cnt[0],
#                     'contact_cnt red het': contact_cnt[3],
#                     'contact_cnt green hom': contact_cnt[1],
#                     'contact_cnt green het': contact_cnt[2],
#                     'Displacement red hom': displacement[0],
#                     'Displacement red het': displacement[3],
#                     'Displacement green hom': displacement[1],
#                     'Displacement green het': displacement[2],
#                     'displacement_var red hom': displacement_var[0],
#                     'displacement_var red het': displacement_var[3],
#                     'displacement_var green hom': displacement_var[1],
#                     'displacement_var green het': displacement_var[2],
#                     'MSD red hom': MSD[0],
#                     'MSD red het': MSD[3],
#                     'MSD green hom': MSD[1],
#                     'MSD green het': MSD[2],
#                     'MSD_var red hom': MSD_var[0],
#                     'MSD_var red het': MSD_var[3],
#                     'MSD_var green hom': MSD_var[1],
#                     'MSD_var green het': MSD_var[2],
#                     })

# df.to_csv(outputFolder+subfolder+"/motility_dataframe.csv")


# Delete zero values and normalise
tmp1 = contact_cnt
tmp2 = distance
tmp3 = displacement
tmp4 = MSD
tmp5 = MSD_var
tmp6 = displacement_var
# contact_cnt=tmp1
# distance=tmp2
# displacement = tmp3
# MSD=tmp4
# Introduce normalised quantities
distance_norm = []
MSD_norm = []
MSD_var_norm = []
displacement_norm = []
displacement_var_norm = []

for idx in range(4):
    cutoff_idx = np.where(contact_cnt[idx,:]>0)   
    distance_norm.append((distance[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])
    MSD_norm.append((MSD[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])
    MSD_var_norm.append((MSD_var[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])
    displacement_norm.append((displacement[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])
    displacement_var_norm.append((displacement_var[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])

# Calculate standard deviation for displacement
displacement_sd = np.zeros((4,cutoff_time))
for idx in range(4):
    displacement_sd[idx,:] = np.sqrt(displacement_var_norm[idx]-displacement_norm[idx])/np.sqrt(contact_cnt[idx,:cutoff_time])

# Calculate standard deviation for MSD
MSD_sd = np.zeros((4,cutoff_time))
for idx in range(4):
    MSD_sd[idx,:] = np.sqrt(MSD_var_norm[idx]-displacement_norm[idx])/np.sqrt(contact_cnt[idx,:cutoff_time])

# Save the plot parameters
fig, ax = plt.subplots(1,2,figsize=(15/2.54, 10/2.54)) # centimeters to inches
ax[0].set_xlabel("Time [min]")
ax[0].set_ylabel(r"Mean squared displacement [$\mu m^2$]")
ax[1].set_xlabel("Time [min]")
ax[1].set_ylabel(r"Total distance travlled [$\mu m$]")

contact_name={0:"EphrinB1 hom",1:"EphB2 hom", 2:"EphB2 het",3:"EphrinB1 het"}
colors = {0:"red",1:"green",2:"green",3:"red"}
markers = {0:'o',1:'o',2:'v',3:'v'}
linestyles={0:'-',1:'-',2:'--',3:'--'}
    

# Fit the Displacement
def f1(t,a):
    return a*t
velocities = np.zeros(4)
velocities_errors = np.zeros(4)
print("Displacements:")
for contact_type in range(4):
    print("Contact type:", contact_name[contact_type])
    measurement_times = range(len(distance_norm[contact_type]))
    X = 2*np.array(measurement_times) # Convert timestep = 2 min
    Y =  3.34*displacement_norm[contact_type] # convert pixel = 3.34 micrometer
    error = 3.34*displacement_sd[contact_type]
    error[0] = np.mean(error)/100 # Can't have zero variance, produces error
    par1, parVar1 = sp.optimize.curve_fit(f1,X, Y, bounds=(0,1e4))#, sigma=error
    velocities[contact_type] = par1[0]
    velocities_errors[contact_type] = np.sqrt(parVar1[0])
    print("Parameters:", par1)
    print("Covariance matrix:", parVar1)
    ax[1].errorbar(X, Y, yerr=error, label=f"{contact_name[contact_type]} ({contact_cnt[contact_type,1]})",color=colors[contact_type],marker=markers[contact_type], ls='none')
    ax[1].plot(X, f1(X, par1[0]), label=f"Fit {contact_name[contact_type]} v = {par1[0]:.3}", color=colors[contact_type],linestyle=linestyles[contact_type])
ax[1].set_title("(a)",  loc='left', fontsize='medium') #fontfamily='sans serif',
D = np.zeros(4)
D_error = np.zeros(4)
tau = np.zeros(4)
tau_error = np.zeros(4)

# # Fit the MSD
# def f2(t,v, tau, D):
#     return (2*D+2*v**2*tau)*t+2*v**2*tau**2*(np.exp(-t/tau)-1)
# print("MSD:")
# for contact_type in range(4):
#     print("Contact type:", contact_type)
#     X = 2*np.array(measurement_times)
#     Y =  4*np.array(MSD_norm[contact_type])
#     error = 4*MSD_sd[contact_type]
#     error[0] = np.mean(error)/100 # Can't have zero variance, produces error
#     # par2, parVar2 = sp.optimize.curve_fit(f2,X, Y,bounds=(0,1e4))
#     # ax[0].plot(X, f2(X, par2[0],par2[1],par2[2]), label=f"Fit {contact_name[contact_type]}  with v={par2[0]:.4}, tau={par2[1]:.4}, D={par2[2]:.4}", color=colors[contact_type],linestyle=linestyles[contact_type])
#     # Fix the velocity
#     f2_fixed = lambda t, tau, D: f2(t,velocities[contact_type], tau, D)
#     par2, parVar2 = sp.optimize.curve_fit(f2_fixed,X, Y,bounds=(0,1e4),sigma=error)
#     ax[0].errorbar(X, Y, yerr=error, label=f"{contact_name[contact_type]} ({contact_cnt[contact_type,1]})",color=colors[contact_type],marker=markers[contact_type], ls='none')
#     ax[0].plot(X, f2_fixed(X, par2[0],par2[1]), label=f"Fit {contact_name[contact_type]}  with v={velocities[contact_type]:.4}, tau={par2[0]:.4}, D={par2[1]:.4}", color=colors[contact_type],linestyle=linestyles[contact_type])
#     print("Parameters:", par2)
#     print("Covariance matrix:", parVar2)

# Fit the MSD without D
def f2(t, tau, v):
    return (2*v**2*tau)*t+2*v**2*tau**2*(np.exp(-t/tau)-1)
print("MSD:")
for contact_type in range(4):
    print("Contact type:", contact_type)
    X = 2*np.array(measurement_times)
    Y =  3.34**2*np.array(MSD_norm[contact_type])
    error = 3.34**2*MSD_sd[contact_type]
    error[0] = np.mean(error)/100 # Can't have zero variance, produces error
    # par2, parVar2 = sp.optimize.curve_fit(f2,X, Y,bounds=(0,1e4))
    # ax[0].plot(X, f2(X, par2[0],par2[1],par2[2]), label=f"Fit {contact_name[contact_type]}  with v={par2[0]:.4}, tau={par2[1]:.4}, D={par2[2]:.4}", color=colors[contact_type],linestyle=linestyles[contact_type])
    # Fix the velocity
    f2_fixed = lambda t, tau: f2(t, tau,velocities[contact_type])
    par2, parVar2 = sp.optimize.curve_fit(f2_fixed,X, Y,bounds=(0,1e4))#,sigma=error
    tau[contact_type] = par2[0]
    tau_error[contact_type] = np.sqrt(parVar2[0])
    # velocities[contact_type] = par2[1]
    # velocities_errors[contact_type] = np.sqrt(parVar2[1,1])
    ax[0].errorbar(X, Y, yerr=error, label=f"{contact_name[contact_type]} ({contact_cnt[contact_type,1]})",color=colors[contact_type],marker=markers[contact_type], ls='none')
    ax[0].plot(X, f2_fixed(X, par2[0]), label=f"Fit {contact_name[contact_type]}  with v={velocities[contact_type]:.3}, tau={par2[0]:.3}", color=colors[contact_type],linestyle=linestyles[contact_type])
    print("Parameters:", par2)
    print("Covariance matrix:", parVar2)

ax[0].legend()
ax[0].set_title("(a)",  loc='left', fontsize='medium') #fontfamily='sans serif',
ax[1].legend()
ax[1].set_title("(b)",  loc='left', fontsize='medium') #fontfamily='sans serif',



# Bar plots for the parameters
fig_par, ax_par = plt.subplots(1,2)
# plt.xticks(rotation = 45)
ax_par[0].set_ylabel(r"Velocity [$\mu m/ min$]")
ax_par[0].bar([0,1,2,3], velocities, yerr=10*velocities_errors, tick_label=["EphrinB1 hom","EphB2 hom","EphB2 het","EphrinB1 het"])
# plt.xticks(rotation = 45)

# fig_tau, ax_tau = plt.subplots(1,1)
ax_par[1].set_ylabel(r"Persistence time [$min$]")
ax_par[1].bar([0,1,2,3], tau, yerr=5*tau_error, tick_label=["EphrinB1 hom","EphB2 hom","EphB2 het","EphrinB1 het"])
# plt.xticks(rotation = 45)
ax_par[0].set_title("(a)",  loc='left', fontsize='medium') #fontfamily='sans serif',
ax_par[1].set_title("(b)",  loc='left', fontsize='medium') #fontfamily='sans serif',
fig_par.autofmt_xdate(rotation=45) # rotate x axis labels
plt.tight_layout() # So labels don't get cut off 

fig.savefig(outputFolder+subfolder+"/motility_errorbars_no_D.png")
fig.savefig(outputFolder+subfolder+"/motility_errorbars_no_D.pdf")
fig_par.savefig(outputFolder+subfolder+"/parameters.pdf")
fig_par.savefig(outputFolder+subfolder+"/parameters.png")
# fig_v.savefig(outputFolder+subfolder+"/velocities_no_D.pdf")
# fig_tau.savefig(outputFolder+subfolder+"/persistence_times_no_D.pdf")

plt.show()
print("Done")