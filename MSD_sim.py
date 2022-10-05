import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
plt.close("all")
"""
Calculates and plots the mean squared displacement for red and green 
particles and compares them to theoretical predictions.
"""

def analytic_MSD(t, Pe, D):
    D_real = D**2 #The D we use in the equations is the square root of the normal radial diffusion constant
    return (4+2*Pe**2/D_real)*t + 2*Pe**2/D_real**2*(np.exp(-t*D_real)-1)


if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/validation/MSD/MSD_parameters"]


for PATH in PATHS:
    green, red, params = load_simulation_data(PATH)

    #Get MSDs
    MSD_green, time = green.MSD()
    MSD_red, time = red.MSD()
    analytic_MSD_green = analytic_MSD(time, params["Pe"], params["greenD"])
    analytic_MSD_red = analytic_MSD(time, params["Pe"], params["redD"])
    
    #plot MSDs
    plt.loglog(time,MSD_green,label="Green (D=%f)"%params["greenD"], color="green")
    plt.loglog(time,MSD_red,label="Red (D=%f)"%params["redD"],color="red")
    plt.loglog(time,analytic_MSD_green, label="Green analytic", color="green", linestyle="--")
    plt.loglog(time,analytic_MSD_red, label="Red analytic", color="red", linestyle="--")

    plt.legend()
    plt.savefig(PATH+"_MSD.png")
    plt.show()
    