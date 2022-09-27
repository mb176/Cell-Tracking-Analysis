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


if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/mixing_index_3/areaFraction_0.4_Pe_40"]
    # PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity/LowDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity/HighDensity",
    #         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl/HighDensityControl"]


#Parameters
neighbourDistance = 1.05

for PATH in PATHS:
    green, red, params = load_simulation_data(PATH)

    # Find last time step
    assert(green.find_t_max()==red.find_t_max()), "Tracks don't finish at the same time"
    t = green.find_t_max()

    #Get mixing index
    mixingIdx, mixing = calculate_mixing_index([green,red], t, neighbourDistance = neighbourDistance)

    print(f"Mixing index: {mixingIdx}")

    #Write mixing index into file
    paramsToWrite = {"Pe":params["Pe"], "areaFraction":params["areaFraction"]}
    strIdx = PATH.rfind("/")
    FILE_PATH = PATH[:strIdx+1]+"mixingIndex"
    write_avr_mixing_index(paramsToWrite, mixingIdx, FILE_PATH)

