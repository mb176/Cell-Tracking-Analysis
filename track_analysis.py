import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from analysis_library import *
from matplotlib import pyplot as plt
from scipy import stats
plt.style.use('seaborn-whitegrid')
plt.close("all")

SOURCE_FOLDER="/home/marius/PhD/CellMotility/tracking_ignacio"
TARGET_FOLDER="/home/marius/PhD/CellMotility/Plots"

names = ['/HighDensitycontrolEphB2/HighDensitycontrolEphB2_greenframes0to211',
        '/HighDensitycontrolEphB2/HighDensitycontrolEphB2_redframes0to211',
        '/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211',
        '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211',
        '/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211',
        '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211'
        ]

########################  Global parameters ##########################
min_length = 0 #Tracks shorter than that are excluded




###################### Singe Cell Type Measurements #######################

# for name in names:
#     exp = experiment(name)
#     exp.read_xml(SOURCE_FOLDER+name+'.xml', min_length)

#     #Plot TAMSD
#     fig, axes = plt.subplots(1,1)
#     exp.plot_TAMSD(tmax, min_length_TAMSD, axes)
#     axes.legend()
#     fig.savefig(TARGET_FOLDER+name+'_TAMSD.svg',format='svg')

#     #Plot Radial distribution function
#     fig, axes = plt.subplots(1,1)
#     n_bins = 150 #Performance bottleneck, since more bins means more randomly generated reference particles
#     for time in [0, 200]:
#         exp.plot_radial_density(axes, time , n_bins, 't = %e'%time)
    
#     axes.legend()
#     fig.savefig(TARGET_FOLDER+name+'_RDF.svg', format='svg)

    # #Plot Tortuosity
    # fig, axes = plt.subplots(1,1)
    # exp.plot_tortuosity(axes, min_length_tortuosuty)
    # fig.savefig(TARGET_FOLDER+name+'_tortuosity.jpg')
    


############################### Mixed Cell Type Measurements ##############################

name_pairs = [['/HighDensitycontrolEphB2/HighDensitycontrolEphB2_greenframes0to211',
        '/HighDensitycontrolEphB2/HighDensitycontrolEphB2_redframes0to211'],
        ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211',
        '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211'],
        ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211',
        '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211']
        ]

for pair in name_pairs:
    #Load green tracks
    green_name = pair[0]
    green = experiment(green_name)
    green.read_xml(SOURCE_FOLDER+green_name+'.xml', min_length)
    
    

    #Load red tracks
    red_name = pair[1]
    red = experiment(red_name)
    red.read_xml(SOURCE_FOLDER+red_name+'.xml', min_length)
    fig, axes = plt.subplots(1,1)

    
    # #Plot both Tortuosities 
    # min_length_tortuosity = 20
    # delta_t = 10
    # n_bins = 25
    # #Dun method
    # fig, axes = plt.subplots(1,1)
    # green.plot_tortuosity(axes, min_length_tortuosity, 'Green', n_bins=n_bins, 
    #                                 method='Dun', delta_t=delta_t)
    # red.plot_tortuosity(axes, min_length_tortuosity, 'Red', 
    #                                 n_bins=n_bins, method='Dun', delta_t=delta_t)
    # axes.set_title(stats.ks_2samp(green.tortuosity,red.tortuosity)) 
    # fig.savefig(TARGET_FOLDER+green_name+'_tortuosity_combined_dun.svg', format='svg')
    # #Normal method
    # fig, axes = plt.subplots(1,1)
    # green.plot_tortuosity(axes, min_length_tortuosity, 'Green', n_bins=n_bins)
    # red.plot_tortuosity(axes, min_length_tortuosity, 'Red', n_bins=n_bins)
    # axes.set_title(stats.ks_2samp(green.tortuosity,red.tortuosity))
    # fig.savefig(TARGET_FOLDER+green_name+'_tortuosity_combined.svg', format='svg')


    # #Plot both radial distributions
    # fig, axes = plt.subplots(1,2)
    # n_bins = 150
    # cutoff_percentage = 70
    # n_reference_points = 10000
    
    # #Plot initial times on the left
    # time = 10
    # green.plot_radial_density(axes[0], time , n_bins, 'Green',cutoff_percentange = cutoff_percentage, 
    #                         n_reference_points = n_reference_points)
    
    # red.plot_radial_density(axes[0], time , n_bins, 'Red',cutoff_percentange = cutoff_percentage, 
    #                         n_reference_points = n_reference_points)
    # axes[0].set_title('t = %e'%time)

    # #Plot final time on the right
    # time = 190
    # green.plot_radial_density(axes[1], time , n_bins, 'Green',cutoff_percentange = cutoff_percentage, 
    #                         n_reference_points = n_reference_points)
    # red.plot_radial_density(axes[1], time , n_bins, 'Red',cutoff_percentange = cutoff_percentage, 
    #                         n_reference_points = n_reference_points)
    # axes[1].set_title('t = %e'%time)
    # fig.savefig(TARGET_FOLDER+green_name+'_RDF_both.svg', format='svg')
    

    # #Plot green-red radial distribution
    n_bins = 150
    cutoff_percentage = 70
    n_reference_points =5000
    fig, axes = plt.subplots(1,2)
    #Early time
    t = 10,
    plot_mixed_particle_radial_density(axes[0], [green,red], t ,n_bins, n_reference_points, no_plot=False, cutoff_percentage=cutoff_percentage)

    #Late time
    t = 190,
    plot_mixed_particle_radial_density(axes[1], [green,red], t ,n_bins, n_reference_points, no_plot=False, cutoff_percentage=cutoff_percentage)
    subfolder = pair[0][0:pair[0][1:].find('/')+1]
    fig.savefig(TARGET_FOLDER+subfolder+'/cross_species_RDF.svg',format='svg')

