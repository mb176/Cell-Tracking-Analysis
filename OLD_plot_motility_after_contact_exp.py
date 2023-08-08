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
contact_radius = 8 # Defines what we mean by contact
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

# Read data frame (the 0th column is the index/ time step)
df = pd.read_csv(outputFolder+subfolder+"/motility_dataframe.csv", index_col=0)

# Use multi-index to distinguish between observables and contact types
df.columns = pd.MultiIndex.from_product([["contact_cnt", "Displacement", "Displacement_var", "MSD", "MSD_var"],["red hom", "red het ", "green hom", "green het"] ], names=['Oberservable', 'Contact type'])

# Filter and normalise results 
df = df.drop(df[df.score < 50].index)

# Using a data frame is actually way less handy than I thought, just filtering 
# zero values with the multiindex is a huge hustle

# Delete zero values and normalise
tmp1 = contact_cnt
tmp2 = distance
tmp3 = displacement
tmp4 = MSD
# contact_cnt=tmp1
# distance=tmp2
# displacement = tmp3
# MSD=tmp4
distance_norm = []
MSD_norm = []
displacement_norm = []
cutoff = 1 # Avoid division by zero
for idx in range(4):
    cutoff_idx = np.where(contact_cnt[idx,:]>cutoff)   
    distance_norm.append((distance[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])
    MSD_norm.append((MSD[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])
    displacement_norm.append((displacement[idx,cutoff_idx]/contact_cnt[idx,cutoff_idx])[0,:cutoff_time])


fig, ax = plt.subplots(1,2) 
ax[0].set_xlabel("Time")
ax[0].set_ylabel("MSD since last contact")
ax[1].set_xlabel("Time")
ax[1].set_ylabel("Displacement since last contact")

contact_name={0:"red-red",1:"green-green", 2:"green with red",3:"red with green"}
colors = {0:"red",1:"green",2:"green",3:"red"}
markers = {0:'o',1:'o',2:'v',3:'v'}
linestyles={0:'-',1:'-',2:'--',3:'--'}
for contact_type in range(4):
    measurement_times = range(len(distance_norm[contact_type]))
    ax[0].scatter(measurement_times, MSD_norm[contact_type], label=f"{contact_name[contact_type]} ({contact_cnt[contact_type,1]})",color=colors[contact_type],marker=markers[contact_type])
    ax[1].scatter(measurement_times, displacement_norm[contact_type], label=f"{contact_name[contact_type]} ({contact_cnt[contact_type,1]})",color=colors[contact_type],marker=markers[contact_type])

# Fit the Displacement
def f1(t,a):
        return a*t
velocities = np.zeros(4)
for contact_type in range(4):
    X = measurement_times
    Y =  displacement_norm[contact_type]
    par1, parVar1 = sp.optimize.curve_fit(f1,X, Y, bounds=(0,1e4))
    velocities[contact_type] = par1[0]
    print("Parameters:", par1)
    print("Covariance matrix:", parVar1)
    ax[1].plot(X, f1(X, par1[0]), label=f"Fit {contact_name[contact_type]} a*x, a = {par1[0]:.4}", color=colors[contact_type],linestyle=linestyles[contact_type])

# Fit the MSD
def f2(t,v, tau, D):
    return (2*D+2*v**2*tau)*t+2*v**2*tau**2*(np.exp(-t/tau)-1)
for contact_type in range(4):
    X = np.array(measurement_times)
    Y =  np.array(MSD_norm[contact_type])
    # par2, parVar2 = sp.optimize.curve_fit(f2,X, Y,bounds=(0,1e4))
    # ax[0].plot(X, f2(X, par2[0],par2[1],par2[2]), label=f"Fit {contact_name[contact_type]}  with v={par2[0]:.4}, tau={par2[1]:.4}, D={par2[2]:.4}", color=colors[contact_type],linestyle=linestyles[contact_type])
    # Fix the velocity
    f2_fixed = lambda t, tau, D: f2(t,velocities[contact_type], tau, D)
    par2, parVar2 = sp.optimize.curve_fit(f2_fixed,X, Y,bounds=(0,1e4))
    ax[0].plot(X, f2_fixed(X, par2[0],par2[1]), label=f"Fit {contact_name[contact_type]}  with v={velocities[contact_type]:.4}, tau={par2[0]:.4}, D={par2[1]:.4}", color=colors[contact_type],linestyle=linestyles[contact_type])
   

    print("Parameters:", par2)
    print("Covariance matrix:", parVar2)

ax[0].legend()
ax[1].legend()
fig.savefig(outputFolder+subfolder+"motility.pdf", format='pdf')

plt.show()
print("Done")