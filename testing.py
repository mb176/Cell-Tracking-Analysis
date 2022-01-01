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
min_length = 0
green = experiment('Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211')
green.read_xml(SOURCE_FOLDER+'/Low Density sorting ephrinB1/Low Density sorting EphB2_green frames 11 to 211'+'.xml', min_length)
red = experiment('Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211')
red.read_xml(SOURCE_FOLDER+'/Low Density sorting ephrinB1/Low Density sorting ephrinB1_red frames 11 to 211'+'.xml', min_length)

fig, axes = plt.subplots(1,2)
t = 10,
n_bins = 60
n_reference_points =3000
label = 'test'
plot_mixed_particle_radial_density(axes[0], [green,red], t ,n_bins, n_reference_points, label, no_plot=False, cutoff_percentage=5)

plt.show()

# green.plot_TAMSD(tmax, min_length_TAMSD, axes[0])


# time = 200
# n_bins = 50
# 
# green.plot_radial_density(axes[0], time , n_bins, 'Green',cutoff_percentange = 10, 
#                             n_reference_points = 1000)#'t = %e'%time
# red.plot_radial_density(axes[0], time , n_bins, 'Red',cutoff_percentange = 10, 
#                            n_reference_points = 100)#'t = %e'%time##
# axes[0].legend()
# plt.show()
# fig.savefig(TARGET_FOLDER+'/test.jpg')
# plt.close('all')