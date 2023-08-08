import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # This avoids crashes on the math cluster from using to many resources
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib
from matplotlib import animation
import copy
import sys
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
matplotlib.rcParams.update({'font.size': 10})
plt.close("all")


"""
Reads all clustering.csv files in the folder and plots the average maximum
cluster size over time for all simulations in one graph (average is over multiple
realisations of the same simulation). 
"""


############################# Parameters #######################################
restrict_parameters = False # If False the selection below is ignored
A_allowed = [0.3, 0.5]
D_allowed = [400, 600, 800]
D_persistent_allowed = [20, 35, 50]
D_SWEEP = True #Set true if parameters are varied over D, D+; Otherwise sweep over Pe, A
D_TAU_SWEEP = False
VARIANT_SWEEP = False # Sweep over different varaints of the model, use the file names as labels
K_A_SWEEP = False
D_A_SWEEP = False
NAME = "avr_cluster_size"
examples = [[10,0.1],[100,0.1],[10,10],[10,50]]
if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output_23_06/CIL_based_demixing/no_cooldown/vary_D_persistentD_0.8/D_600_tau_0.001"


# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")+1] # rindex finds last occurence of a character 


# Initialise graphic
fig_all, ax_all = plt.subplots(1,1)
fig_colors, ax_colors = plt.subplots(1,2)
ax_green = ax_colors[0]
ax_red = ax_colors[1]
fig_all_log, ax_all_log = plt.subplots(1,1)
fig_colors_log, ax_colors_log = plt.subplots(1,2)
ax_green_log = ax_colors_log[0]
ax_red_log = ax_colors_log[1]

# Plotting style
lineStyles = ["solid", (0,(10,1)), (0,(8,3)), (0,(7,4)), (0,(5,5)), (0,(3,3)), (0,(1,1)), (0,(1,2)), (0,(1,1,3,1))]
lineColors = ["black", "blue", "orange", "green", "red", "purple", "brown", "gray", "olive", "cyan", "pink"]
styleIdx = 0
colorIdx = 0
A_dic = {}
Pe_dic = {}
D_dic = {}
tau_dic = {}
persistentD_dic = {}
k_dic = {}

# Prepare phasediagramm of the mixing index in the last timestep for all particles
fig_phasediagramm, axes = plt.subplot_mosaic([['a)', 'b)', 'c)'], ['a)','d)', 'e)']],
                              layout='constrained', 
                              gridspec_kw={'width_ratios': [2,1,1]},
                              figsize = (16/2.54,8/2.54))
ax2 = axes['a)']
# fig_parameter_plot, ax3 = plt.subplots(1,1,figsize=(10/2.54, 10/2.54))

areaFractions = []
Pes = []
Ds = []
persistentDs = []
taus = []
ks = []
values = []
valuesRed = []
valuesGreen = []




for fileName in sorted(os.listdir(FOLDER)): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+fileName[:-11]


    if paramFile is not None and paramFile[-10:]!="velocities":
        INPUT_FILE = paramFile+"_clustering.csv"
        # Find parameters
        exp = experiment(paramFile)
        params = exp.read_parameter_file(paramFile)
        A = params["areaFraction"]
        Pe = params["Pe"]
        nParticles = params["nParticles"]
        nGreenParticles = params["nGreenParticles"]
        nRedParticles = params["nRedParticles"]
        D = params["greenD"]
        persistentD = params["greenPersistentD"]
        tau = params["tau"]
        k = params["k"]


        if (k>= 1000) or restrict_parameters==False:
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

            # Read cluster sizes from file
            clusterSize, clusterSizeGreen, clusterSizeRed, measurementTimes = read_clustering_file(INPUT_FILE)
            # # Plot  cluster size histogram for each experiment
            # fig_histogram, ax1 = plt.subplots(1,1)
            # tIdx = -1
            # plot_cluster_histogram(ax1, clusterSizeFrequencies[tIdx], label=f"A={A}, Pe={Pe}" ,color="green")        
            # ax1.legend()
            # fig_histogram.savefig(paramFile+"cluster_histogram.png", format='png',dpi=300)

            # Average maximum cluster size over all realisations
            returnAverage = True
            avrMaxClusterSize = get_avrMaxClusterSize(measurementTimes, clusterSize, nParticles, average=returnAverage)
            avrMaxClusterSizeGreen = get_avrMaxClusterSize(measurementTimes, clusterSizeGreen, nGreenParticles, average=returnAverage)
            avrMaxClusterSizeRed = get_avrMaxClusterSize(measurementTimes, clusterSizeRed, nRedParticles, average=returnAverage)

            # maxClusterSize = [max(clusterSize[tIdx]) for tIdx in range(len(measurementTimes))]
            # maxClusterSizeGreen = [max(clusterSizeGreen[tIdx]) for tIdx in range(len(measurementTimes))]
            # maxClusterSizeRed = [max(clusterSizeRed[tIdx]) for tIdx in range(len(measurementTimes))]

            if D_SWEEP:
                ax_all.plot(measurementTimes, avrMaxClusterSize, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
                ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
                ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
                ax_all_log.loglog(measurementTimes, avrMaxClusterSize, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
                ax_green_log.loglog(measurementTimes, avrMaxClusterSizeGreen, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
                ax_red_log.loglog(measurementTimes, avrMaxClusterSizeRed, label=f"D={D}, D+={persistentD}", color=D_dic[D], linestyle=persistentD_dic[persistentD])
            elif VARIANT_SWEEP:
                ax_all.plot(measurementTimes, avrMaxClusterSize, label=fileName[:-13], color=lineColors[colorIdx])
                ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=fileName[:-13], color=lineColors[colorIdx])
                ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=fileName[:-13], color=lineColors[colorIdx])
                ax_all_log.loglog(measurementTimes, avrMaxClusterSize, label=fileName[:-13], color=lineColors[colorIdx])
                ax_green_log.loglog(measurementTimes, avrMaxClusterSizeGreen, label=fileName[:-13], color=lineColors[colorIdx])
                ax_red_log.loglog(measurementTimes, avrMaxClusterSizeRed, label=fileName[:-13], color=lineColors[colorIdx])
                colorIdx+=1
            elif D_TAU_SWEEP:
                ax_all.plot(measurementTimes, avrMaxClusterSize, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
                ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
                ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
                ax_all_log.loglog(measurementTimes, avrMaxClusterSize, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
                ax_green_log.loglog(measurementTimes, avrMaxClusterSizeGreen, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
                ax_red_log.loglog(measurementTimes, avrMaxClusterSizeRed, label=f"D={D}, tau={tau}", color=D_dic[D], linestyle=tau_dic[tau])
            elif D_A_SWEEP:
                ax_all.plot(measurementTimes, avrMaxClusterSize, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
                ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
                ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
                ax_all_log.loglog(measurementTimes, avrMaxClusterSize, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
                ax_green_log.loglog(measurementTimes, avrMaxClusterSizeGreen, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
                ax_red_log.loglog(measurementTimes, avrMaxClusterSizeRed, label=f"D={D}, A={A}", color=D_dic[D], linestyle=A_dic[A])
            elif K_A_SWEEP:
                ax_all.plot(measurementTimes, avrMaxClusterSize, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
                ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
                ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
                ax_all_log.loglog(measurementTimes, avrMaxClusterSize, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
                ax_green_log.loglog(measurementTimes, avrMaxClusterSizeGreen, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
                ax_red_log.loglog(measurementTimes, avrMaxClusterSizeRed, label=f"k={k}, A={A}", color=k_dic[k], linestyle=A_dic[A])
            else:
                ax_all.plot(measurementTimes, avrMaxClusterSize, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                ax_green.plot(measurementTimes, avrMaxClusterSizeGreen, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                ax_red.plot(measurementTimes, avrMaxClusterSizeRed, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                ax_all_log.loglog(measurementTimes, avrMaxClusterSize, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                ax_green_log.loglog(measurementTimes, avrMaxClusterSizeGreen, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])
                ax_red_log.loglog(measurementTimes, avrMaxClusterSizeRed, label=f"A={A}, Pe={Pe}", color=Pe_dic[Pe], linestyle=A_dic[A])

            # Save final mixing index for phase diagram
            tIdx = -1
            areaFractions.append(A)
            Pes.append(Pe)
            Ds.append(D)
            taus.append(tau)
            persistentDs.append(persistentD)
            ks.append(k)
            values.append(avrMaxClusterSize[-1])
            valuesRed.append(avrMaxClusterSizeRed[-1])
            valuesGreen.append(avrMaxClusterSizeGreen[-1])

# Label and save the plots
ax_all.legend(); ax_all.set_xlabel("Time"); ax_all.set_ylabel("Average cluster size"); ax_all.set_title("All clusters")
ax_green.legend(); ax_green.set_xlabel("Time"); ax_green.set_ylabel("Average cluster size"); ax_green.set_title("Green clusters")
ax_red.legend(); ax_red.set_xlabel("Time"); ax_red.set_ylabel("Average cluster size"); ax_red.set_title("Red clusters")
ax_all_log.legend(); ax_all_log.set_xlabel("Time"); ax_all_log.set_ylabel("Average cluster size"); ax_all_log.set_title("All clusters")
ax_green_log.legend(); ax_green_log.set_xlabel("Time"); ax_green_log.set_ylabel("Average cluster size"); ax_green_log.set_title("Green clusters")
ax_red_log.legend(); ax_red_log.set_xlabel("Time"); ax_red_log.set_ylabel("Average cluster size"); ax_red_log.set_title("Red clusters")

# Make room for the legend on the left
ax_all.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_green.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_red.set_xlim(-0.5*measurementTimes[-1], measurementTimes[-1]) 
ax_all_log.set_xlim(measurementTimes[1], 1.5*measurementTimes[-1]) 
ax_red_log.set_xlim(measurementTimes[1], 1.5*measurementTimes[-1]) 
ax_green_log.set_xlim(measurementTimes[1], 1.5*measurementTimes[-1]) 

# fig_green.savefig(FOLDER+NAME+"_green.pdf", format='pdf')
# fig_red.savefig(FOLDER+NAME+"_red.pdf", format='pdf')
# fig_all_log.savefig(FOLDER+NAME+"_log.pdf", format='pdf')
# fig_green_log.savefig(FOLDER+NAME+"_green_log.pdf", format='pdf')
# fig_red_log.savefig(FOLDER+NAME+"_red_log.pdf", format='pdf')
#Fix color map:

# norm=matplotlib.colors.Normalize(min(min(valuesGreen),min(valuesRed)), max(max(valuesGreen),max(valuesRed)))
# Map to standard definition of rotational diffusion coefficient
Ds = np.array(Ds)**2
persistentDs = np.array(persistentDs)**2
val = np.array(valuesGreen)-np.array(valuesRed)

norm=matplotlib.colors.CenteredNorm() # Center 0 on white min(val), max(val)
if D_SWEEP:
    # for D in Ds:
    #     X = [persistentDs[idx] for idx in range(len(val)) if Ds[idx]==D]
    #     Y = [val[idx] for idx in range(len(val)) if Ds[idx]==D]
    #     ax3.plot(X, Y, label=f"D={D}")
    # Get divergent color map "PuOr" (two colors)
    im1 = ax2.scatter(Ds, persistentDs, s=200, c=np.array(valuesGreen)-np.array(valuesRed), norm=norm, cmap='PuOr') #cmap='Blues'
    # im2 = ax2[1].scatter(Ds, persistentDs, s=200, c=valuesRed, cmap='Blues', norm=norm)
    ax2.set_xlabel(r"$D$")
    ax2.set_ylabel(r"$D^+$")
    # ax2[1].set_xlabel(r"$D$")
    # ax2[1].set_ylabel(r"$D_+$")
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
elif D_TAU_SWEEP:
    im1 = ax2.scatter(Ds, taus, s=200, c=np.array(valuesGreen)-np.array(valuesRed), cmap='PuOr', norm=norm)
    # im2 = ax2[1].scatter(Ds, taus, s=200, c=valuesRed, cmap='Blues', norm=norm)
    ax2.set_xlabel(r"$D$")
    ax2.set_ylabel(r"$\tau$")
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    # ax2[1].set_xlabel(r"$D$")
    # ax2[1].set_ylabel(r"$\tau$")
elif D_A_SWEEP:
    im1 = ax2.scatter(Ds, areaFractions, s=200, c=np.array(valuesGreen)-np.array(valuesRed), cmap='PuOr', norm=norm)
    # im2 = ax2[1].scatter(Ds, areaFractions, s=200, c=valuesRed, cmap='Blues', norm=norm)
    ax2.set_xlabel(r"$D$")
    ax2.set_ylabel(r"Area Fraction")
    # ax2[1].set_xlabel(r"$D$")
    # ax2[1].set_ylabel(r"Area Fraction")
elif K_A_SWEEP:
    im1 = ax2.scatter(ks, areaFractions, s=200, c=np.array(valuesGreen)-np.array(valuesRed), cmap='PuOr', norm=norm)
    # im2 = ax2[1].scatter(ks, areaFractions, s=200, c=valuesRed, cmap='Blues', norm=norm)
    ax2.set_xscale("log")
    # ax2[1].set_xscale("log")
    ax2.set_xlabel(r"$k$")
    ax2.set_ylabel(r"Area Fraction")
    # ax2[1].set_xlabel(r"$k$")
    # ax2[1].set_ylabel(r"Area Fraction")
else:
    im1 = ax2.scatter(areaFractions, Pes, s=200, c=np.array(valuesGreen)-np.array(valuesRed), cmap='PuOr', norm=norm)
    # im2 = ax2[1].scatter(areaFractions, Pes, s=200, c=valuesRed, cmap='Blues', norm=norm)
    ax2.set_xlabel("Area Fraction")
    ax2.set_ylabel("Peclet")
    # ax2[1].set_xlabel("Area Fraction")
    # ax2[1].set_ylabel("Peclet")


# Colorbar for phase diagram
# fig_phasediagramm.subplots_adjust(right=0.8)
# cbar_ax = fig_phasediagramm.add_axes([0.85, 0.15, 0.05, 0.7])
# fig_phasediagramm.colorbar(im1, cax=cbar_ax)
ax2.set_title(r"$<N_{green}>-<N_{red}>$",fontsize='10')

# Add the Examples
ax2.set_title("(a)",  loc='left', fontsize='10') #fontfamily='sans serif',
axes_labels = ["b)","c)","d)","e)"]

for idx in range(4):
    ax = axes[axes_labels[idx]]
    PATH = FOLDER +f"D_{examples[idx][0]}_persistentD_{examples[idx][1]}"
    final_snapshot(ax,PATH, velocityArrows=False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title("("+axes_labels[idx],  loc='left', fontsize='10') #fontfamily='sans serif',
    # ax.set_title(f"$D$={examples[idx][0]}, $D^+$={examples[idx][1]}" ,  loc='right', fontsize='8') #fontfamily='sans serif',
    # ax.set_title(f"A={examples[idx][0]}, D={examples[idx][0]}",  loc='center', fontsize='medium') #fontfamily='sans serif',

# ax2.set_aspect(0.8) # make room for colorbar
fig_phasediagramm.colorbar(im1,ax=ax2)#,shrink=0.6
# plt.tight_layout()
fig_phasediagramm.savefig(FOLDER+"clustering_examples.png",dpi=400)

#Save figures
fig_all.savefig(FOLDER+NAME+".png")
fig_colors.savefig(FOLDER+NAME+"_by_type.png")
fig_colors_log.savefig(FOLDER+NAME+"_by_typ_log.png")
# fig_phasediagramm.savefig(FOLDER+NAME+"_phasediagram.png")