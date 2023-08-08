import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from analysis_library import *
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
plt.style.use('seaborn-whitegrid')
import matplotlib
matplotlib.rcParams.update({'font.size': 10})
plt.rcParams["font.family"] = "sans serif"
plt.close("all")

SOURCE_FOLDER="/home/marius/PhD/CellMotility/tracking_23_01" #tracking_ignacio_2022"   #  
TARGET_FOLDER="/home/marius/PhD/CellMotility/Plots/Plots_2023_01" #Plots_Ignacio2022"   #

########################  Global parameters ##########################
min_length = 0 #Tracks shorter than that are excluded



# name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
#         '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
#         ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
#         '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks'],
#         ['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
#         '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']
#         ]


# High density
name_pairs = [['/HighDensitycontrolEphB2/High Density control EphB2_green frames 0 to 211_Tracks',
        '/HighDensitycontrolEphB2/High Density control EphB2_red frames 0 to 211_Tracks'],
        ['/High Density sorting EphB2/High Density sorting EphB2_green frames 0 to 211_Tracks',
        '/High Density sorting EphB2/High Density sorting ephrinB1_red 0 to 211_Tracks']]

# [['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211',
#                 '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211']]

n_bins = 200 
cutoff_percentage = 80
n_reference_points = 2000
times = [10,190]#np.linspace(0,190,20) 
fig, ax = plt.subplots(2,2,figsize=(15/2.54, 15/2.54),sharey=True) # centimeters to inches

title=["Control","Sorting"]

for pairIdx, pair in enumerate(name_pairs):
    
    #Load green tracks
    green_name = pair[0][:pair[0].find('f')]
    green = experiment(green_name)
    green.read_xml(SOURCE_FOLDER+pair[0]+'.xml', min_length)
    
    

    #Load red tracks
    red_name = pair[1][:pair[1].find('f')]
    red = experiment(red_name)
    red.read_xml(SOURCE_FOLDER+pair[1]+'.xml', min_length)
    

    
    timeIdx = 0
    for time in times:
        axes = ax[pairIdx, timeIdx]
        timeIdx += 1
        #Green particles
        green.plot_radial_density(axes, time , n_bins, 'EphB2','g',cutoff_percentange = cutoff_percentage, 
                                n_reference_points = n_reference_points)
        #red particles
        red.plot_radial_density(axes, time , n_bins, 'EphrinB1','r',cutoff_percentange = cutoff_percentage, 
                                n_reference_points = n_reference_points)
        #red-green crosscorrelation
        plot_mixed_particle_radial_density(axes, [green,red], time ,n_bins, n_reference_points, 
                                            cutoff_percentage=cutoff_percentage)
        bin_size = np.sqrt(green.x_max**2+green.y_max**2)/n_bins
        axes.set_title(title[pairIdx]+', t = %d min'%(2*time))

        
        # subfolder = green_name[:green_name.find('/',1,len(green_name))]
    #     
plt.tight_layout()
fig.savefig(TARGET_FOLDER+'/RDF/'+'High_density_RDF_d_%f.png'%(bin_size), format='png', dpi=300)




# Low Desnity
name_pairs = [['/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211_Tracks',
                '/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211_Tracks']]


fig, ax = plt.subplots(1,2,figsize=(15/2.54, 9/2.54),sharey=True) # centimeters to inches

for pairIdx, pair in enumerate(name_pairs):
    
    #Load green tracks
    green_name = pair[0][:pair[0].find('f')]
    green = experiment(green_name)
    green.read_xml(SOURCE_FOLDER+pair[0]+'.xml', min_length)
    
    

    #Load red tracks
    red_name = pair[1][:pair[1].find('f')]
    red = experiment(red_name)
    red.read_xml(SOURCE_FOLDER+pair[1]+'.xml', min_length)
    
    
    
    timeIdx = 0
    for time in times:
        axes = ax[timeIdx]
        timeIdx += 1
        #Green particles
        green.plot_radial_density(axes, time , n_bins, 'EphB2','g',cutoff_percentange = cutoff_percentage, 
                                n_reference_points = n_reference_points)
        #red particles
        red.plot_radial_density(axes, time , n_bins, 'EphrinB1','r',cutoff_percentange = cutoff_percentage, 
                                n_reference_points = n_reference_points)
        #red-green crosscorrelation
        plot_mixed_particle_radial_density(axes, [green,red], time ,n_bins, n_reference_points, 
                                            cutoff_percentage=cutoff_percentage)
        bin_size = np.sqrt(green.x_max**2+green.y_max**2)/n_bins
        axes.set_title('t = %d min'%(2*time))

        
        # subfolder = green_name[:green_name.find('/',1,len(green_name))]
    #     
plt.tight_layout()
fig.savefig(TARGET_FOLDER+'/RDF/'+'low_density_RDF_d_%f.png'%(bin_size), format='png')

plt.show()