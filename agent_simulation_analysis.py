import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
plt.close("all")


#Parameter file given externally:
PATHS = [sys.argv[1]]

#Give parameter file manually
# PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters"]
# PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity/LowDensity",
#         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity/HighDensity",
#         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl/HighDensityControl"]



for PATH in PATHS:
    #Load the data
    green = experiment(PATH)
    green.read_csv(PATH+"_tracks.csv")
    params = green.read_parameter_file(PATH)


    # Global parameters
    min_length = 0 #Tracks shorter than that are excluded
    nGreenParticles = int(params["nGreenParticles"])
    nRedParticles = int(params["nRedParticles"])

    #split up red and green particles (the original experiment will store green tracks)
    red = copy.deepcopy(green)
    green.tracks = green.tracks[:nGreenParticles]
    green.color = green.color[:nGreenParticles]
    red.tracks = red.tracks[nGreenParticles:nGreenParticles+nRedParticles]
    red.color = red.color[nGreenParticles:nGreenParticles+nRedParticles]

    # #Plot both Tortuosities 
    # min_length_tortuosity = 10
    # delta_t = 5 #time interval for Dun method
    # n_bins = 30

    # #Dun method 
    # fig, axes = plt.subplots(1,1)
    # green.plot_tortuosity(axes, min_length_tortuosity, 'Green cells','g', n_bins=n_bins, 
    #                                 method='Dun', delta_t=delta_t)
    # red.plot_tortuosity(axes, min_length_tortuosity, 'Red cells', 'r',
    #                                 n_bins=n_bins, method='Dun', delta_t=delta_t)
    # axes.set_title(stats.ks_2samp(green.tortuosity,red.tortuosity)) 
    # fig.savefig(PATH+'_tortuosity_combined_dun_t_%d_min_t_%d.png'%(delta_t, min_length_tortuosity), format='png',dpi=200)

    # #Normal method
    # fig, axes = plt.subplots(1,1)
    # green.plot_tortuosity(axes, min_length_tortuosity, 'Green cells','g', n_bins=n_bins)
    # red.plot_tortuosity(axes, min_length_tortuosity, 'Red cells','r', n_bins=n_bins)
    # axes.set_title(stats.ks_2samp(green.tortuosity,red.tortuosity))
    # fig.savefig(TARGET_FOLDER+green_name+'_tortuosity_combined.svg', format='svg')

    # #Plot RDF Slideshow
    # n_bins = 300 
    # cutoff_percentage = 1
    # n_reference_points = 2000
    # times = np.linspace(0,90,10)
    # for time in times:
    #     fig, axes = plt.subplots(1,1)
    #     #Green particles
    #     green.plot_radial_density(axes, time , n_bins, 'Green cells','g',cutoff_percentange = cutoff_percentage, 
    #                             n_reference_points = n_reference_points)
    #     #red particles
    #     red.plot_radial_density(axes, time , n_bins, 'Red cells','r',cutoff_percentange = cutoff_percentage, 
    #                             n_reference_points = n_reference_points)
    #     # #red-green crosscorrelation
    #     plot_mixed_particle_radial_density(axes, [green,red], time ,n_bins, n_reference_points, 
    #                                         cutoff_percentage=cutoff_percentage)
    #     bin_size = np.sqrt(green.x_max**2+green.y_max**2)/n_bins
    #     axes.set_title('t = %d, bin size = %f'%(time,bin_size))
    #     fig.savefig(PATH+'RDF_t_%i_test.png'%time, format='png',dpi=300)
    

    plt.close("all")

