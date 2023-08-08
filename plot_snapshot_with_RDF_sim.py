import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # This avoids crashes on the math cluster from using to many resources
import matplotlib.pyplot as plt 
import numpy as np
import csv
import os
import matplotlib.animation as animation
from celluloid import Camera
from matplotlib.collections import PatchCollection
import sys
sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
import matplotlib
matplotlib.rcParams.update({'font.size': 10})
plt.close("all")

"""
Plots and saves a picture of the final arrangement of particles.
"""

if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output_23_06/persistence_based_demixing/vary_tau_D/D_600_tau_0.001"]
    # PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity/LowDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity/HighDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl/HighDensityControl"]

for PATH in PATHS:
    fig, ax = plt.subplots(1,2, figsize=(15/2.54, 10/2.54), gridspec_kw={'width_ratios': [2, 1.1]})#

    # Final Snapshot
    final_snapshot(ax[0],PATH, velocityArrows=False)

    
    # RDF
    
    try: # Single realisation
        green, red, params = load_simulation_data(PATH)
    except: # Multiple realisations
        green, red, params = load_simulation_data(PATH,trackIdx="_1")
    
    n_bins = 200 
    cutoff_percentage = 80
    n_reference_points = 1000
    t_final = green.tracks[0][-1][0]
    # x = np.linspace(0,10, 50)
    # ax[1].plot(x, np.sin(x))
    #Green particles
    green.plot_radial_density(ax[1], t_final , n_bins, 'Green','g',cutoff_percentange = cutoff_percentage, 
                            n_reference_points = n_reference_points)
    #red particles
    red.plot_radial_density(ax[1], t_final , n_bins, 'Red','r',cutoff_percentange = cutoff_percentage, 
                            n_reference_points = n_reference_points)
    # #red-green crosscorrelation
    plot_mixed_particle_radial_density(ax[1], [green,red], t_final ,n_bins, n_reference_points, 
                                        cutoff_percentage=cutoff_percentage, label="Green-red")
    # ax_inset.set_ylabel("RDF")
    # ax_inset.set_xlabel("")
    # ax_inset.get_legend().remove()
    
ax[0].set_title("(a)",  loc='left', fontsize='medium') #fontfamily='sans serif',
ax[1].set_title("(b)",  loc='left', fontsize='medium') #fontfamily='sans serif',
# plt.show()
fig.savefig(PATH+"_final_frame_RDF.png",dpi=500)