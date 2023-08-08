import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # This avoids crashes on the math cluster from using to many resources
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import matplotlib
import copy
import sys
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
matplotlib.rcParams.update({'font.size': 10})
plt.close("all")

"""
Reads all mixing.csv files in the folder and plots the average mixing 
index over time for all simulations in one graph (average is over multiple
realisations of the same simulation), also adds snapshots for four selected 
parameter values
"""

############################# Parameters #######################################
# By default the file sweeps over all combinations of A and Pe in the folder
D_SWEEP = False #Set true if parameters are varied over D, D+ instead
D_TAU_SWEEP = False
VARIANT_SWEEP = False # Sweep over different varaints of the model, use the file names as labels
K_A_SWEEP = False
D_A_SWEEP = True
A_allowed = [0.8]
restrict_parameters = False # If False A_allowed is ignored
INPUT_FILE_NAME = "_mixing_d_1.1.csv" #"_mixing_d_1.20978631683.csv" # #"_mixing.csv"" 
OUTPUT_FILE_NAME = "mixing_over_time"
if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output_23_06/CIL_based_demixing/no_cooldown_no_persistence/vary_D_A/A_0.2_D_1"
    
# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 

# Plotting style
lineStyles = ["solid", (0,(6,1)), (0,(3,3)), (0,(6,1,1)), (0,(6,1,1,1)), (0,(6,1,1,1,1,1)), (0,(6,1,1,1,1,1,1,1))]#["solid", (0,(10,1)), (0,(8,3)), (0,(7,4)), (0,(5,5)), (0,(3,3)), (0,(1,1)), (0,(1,2)), (0,(1,1,3,1))]
lineColors = ["black", "blue", "orange", "green", "red", "purple", "brown", "gray", "olive", "cyan", "pink", "magenta"]
styleIdx = 0
colorIdx = 0
A_dic = {}
Pe_dic = {}
D_dic = {}
tau_dic = {}
persistentD_dic = {}
k_dic = {}


# Initialise figures
fig_time_dependence, ax_all = plt.subplots(1,1) 
fig_log, ax_log = plt.subplots(1,1)
# fig_time_dependence_green, ax_green= plt.subplots(1,1)
# fig_time_dependence_red, ax_red = plt.subplots(1,1)

# Prepare phasediagramm of the mixing index in the last timestep for all particles
fig_phasediagramm, axes = plt.subplot_mosaic([['a)', 'b)', 'c)'], ['a)','d)', 'e)']],
                              layout='constrained', 
                              gridspec_kw={'width_ratios': [2,1,1]},
                              figsize = (15/2.54,7.5/2.54))
ax2 = axes['a)']
areaFractions = []
Pes = []
Ds = []
persistentDs = []
taus = []
ks = []
values = []


for fileName in sorted(os.listdir(FOLDER)): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None                
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+fileName[:-11]

    if paramFile is not None and paramFile[-10:]!="velocities":
        
        INPUT_FILE = paramFile+INPUT_FILE_NAME
        if not os.path.exists(INPUT_FILE): # Check if mixing file has been created
            continue
        # Find parameters
        exp = experiment(paramFile)
        params = exp.read_parameter_file(paramFile)
        A = params["areaFraction"]
        Pe = params["Pe"]
        D = params["greenD"]
        persistentD = params["greenPersistentD"]
        tau = params["tau"]
        k = params["k"]

        if (k in [1000, 5000]) or restrict_parameters==False:
            # Update linestyle and line color
            if D_SWEEP:
                if D not in D_dic:
                    D_dic[D]= lineColors[colorIdx]
                    colorIdx += 1
                if persistentD not in persistentD_dic:
                    persistentD_dic[persistentD]=lineStyles[styleIdx]
                    styleIdx += 1
            elif D_TAU_SWEEP:
                if D not in D_dic:
                    D_dic[D]= lineColors[colorIdx]
                    colorIdx += 1
                if tau not in tau_dic:
                    tau_dic[tau]=lineStyles[styleIdx]
                    styleIdx += 1
            elif D_A_SWEEP:
                if D not in D_dic:
                    D_dic[D]= lineColors[colorIdx]
                    colorIdx += 1
                if A not in A_dic:
                    A_dic[A]=lineStyles[styleIdx]
                    styleIdx += 1
            elif K_A_SWEEP:
                if k not in k_dic:
                    k_dic[k]= lineColors[colorIdx]
                    colorIdx += 1
                if A not in A_dic:
                    A_dic[A]=lineStyles[styleIdx]
                    styleIdx += 1
            else:
                if A not in A_dic:
                    A_dic[A]=lineStyles[styleIdx]
                    styleIdx += 1
                if Pe not in Pe_dic:
                    Pe_dic[Pe]= lineColors[colorIdx]
                    colorIdx += 1
            
            

            # Read mixing index from file
            mixingIdx, mixingIdxGreen, mixingIdxRed, measurementTimes = read_mixing_index_file(INPUT_FILE)
            
            # # Plot mixing index over time
            if D_SWEEP:
                ax_all.plot(measurementTimes, mixingIdx, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
                ax_log.loglog(measurementTimes, mixingIdx, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
            elif VARIANT_SWEEP:
                ax_all.plot(measurementTimes, mixingIdx, label=fileName[:-13])#, color=lineColors[colorIdx]
                ax_log.loglog(measurementTimes, mixingIdx, label=fileName[:-13])#, color=lineColors[colorIdx]
                # colorIdx+=1
            elif D_TAU_SWEEP:
                ax_all.plot(measurementTimes, mixingIdx, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
                ax_log.loglog(measurementTimes, mixingIdx, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
            elif D_A_SWEEP:
                ax_all.plot(measurementTimes, mixingIdx, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
                ax_log.loglog(measurementTimes, mixingIdx, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
            elif K_A_SWEEP:
                ax_all.plot(measurementTimes, mixingIdx, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
                ax_log.loglog(measurementTimes, mixingIdx, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
            else:
                ax_all.plot(measurementTimes, mixingIdx, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                ax_log.loglog(measurementTimes, mixingIdx, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                

            # Save final mixing index for phase diagram
            tIdx = -1
            areaFractions.append(A)
            Pes.append(Pe)
            Ds.append(D)
            taus.append(tau)
            persistentDs.append(persistentD)
            ks.append(k)
            values.append(mixingIdx[-1])
        
ax_all.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) # Make room for legend on the left
ax_all.set_xlabel("Time")
ax_all.set_ylabel("Demixing index")
ax_all.legend()
# fig_time_dependence.savefig(FOLDER+OUTPUT_FILE_NAME+".png", format='png')

ax_log.set_xlim(measurementTimes[1], 5*measurementTimes[-1]) 
ax_log.set_xlabel("Time")
ax_log.set_ylabel("Demixing index")
ax_log.legend()
# fig_log.savefig(FOLDER+OUTPUT_FILE_NAME+"_log.png")

# ax_green.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) # Make room for legend on the left
# ax_green.set_xlabel("Time")
# ax_green.set_ylabel("Demixing index")
# ax_green.legend()
# fig_time_dependence_green.savefig(FOLDER+NAME+"_green.pdf", format='pdf')

# ax_red.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) # Make room for legend on the left
# ax_red.set_xlabel("Time")
# ax_red.set_ylabel("Demixing index")
# ax_red.legend()
# fig_time_dependence_red.savefig(FOLDER+NAME+"_red.pdf", format='pdf')
Ds = np.array(Ds)**2
persistentDs = np.array(persistentDs)**2
if D_SWEEP:
    im = ax2.scatter(Ds, persistentDs, s=200, c=values, cmap='Blues')
    ax2.set_xlabel(r"$D$")
    ax2.set_ylabel(r"$D^+$")
elif D_TAU_SWEEP:
    im = ax2.scatter(Ds, taus, s=200, c=values, cmap='Blues')
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$D$")
    ax2.set_ylabel(r"$\tau_p$")
elif D_A_SWEEP:
    im = ax2.scatter(Ds, areaFractions, s=200, c=values, cmap='Blues')
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$D$")
    ax2.set_ylabel(r"Packing Fraction")

elif K_A_SWEEP:
    im = ax2.scatter(ks, areaFractions, s=200, c=values, cmap='Blues')
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$k$")
    ax2.set_ylabel(r"Packing Fraction")
else:
    im = ax2.scatter(areaFractions, Pes, s=200, c=values, cmap='Blues')
    ax2.set_yscale("log")
    ax2.set_xlabel("Packing Fraction")
    ax2.set_ylabel("Peclet")

    # # inset axes....
    # axins = ax2.inset_axes([1.1, 0, 0.45, 0.45])
    # x = np.linspace(0,5,20)
    # axins.plot(x,np.sin(x))
    # # # subregion of the original image
    # # x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9
    # # axins.set_xlim(x1, x2)
    # # axins.set_ylim(y1, y2)
    # axins.set_xticklabels([])
    # axins.set_yticklabels([])
    # plt.tight_layout()

    # ax2.indicate_inset_zoom(axins, edgecolor="black")


# Add the Examples
ax2.set_title("(a)",  loc='left', fontsize='medium') #fontfamily='sans serif',
axes_labels = ["b)","c)","d)","e)"]
examples = [[0.8,1],[0.8,50],[0.4,1],[0.8,1000]]
for idx in range(4):
    ax = axes[axes_labels[idx]]
    PATH = FOLDER +f"A_{examples[idx][0]}_D_{examples[idx][1]}"
    final_snapshot(ax,PATH, velocityArrows=False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title("("+axes_labels[idx],  loc='left', fontsize='10') #fontfamily='sans serif',f" A={examples[idx][0]}, D={examples[idx][0]}" 
    ax.set_title(f"A={examples[idx][0]}, D={examples[idx][1]**2:.1g}",  loc='right', fontsize='8') #fontfamily='sans serif',

# ax2.set_aspect(0.8) # make room for colorbar
fig_phasediagramm.colorbar(im,ax=ax2)#,shrink=0.6
# plt.tight_layout()
fig_phasediagramm.savefig(FOLDER+"mixing_examples.png",dpi=400)

# ax_all.plot(measurementTimes, mixingIdx, label=f"Differential persistence", color=lineColors[colorIdx], linestyle=lineStyles[styleIdx])
# ax_log.loglog(measurementTimes, mixingIdx, label=f"Differential persistence", color=lineColors[colorIdx], linestyle=lineStyles[styleIdx])
# styleIdx+=1; colorIdx+=1