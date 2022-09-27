import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
import os
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *
plt.style.use('seaborn-whitegrid')
plt.close("all")


if(len(sys.argv)==2): #Parameter file given externally:
    PATH = sys.argv[1]
else: #Give parameter file manually
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/new_output/test/test_file_iteration/new_test_parameters"

# Folder containing the PATH
FOLDER = PATH[0:PATH.rindex("/")] # rindex finds last occurence of a character 

# Parameter
NEIGHBOUR_DISTANCE=None

for fileName in os.listdir(FOLDER): #Iterate over all files in the folder
    # Find all parameter files
    paramFile = None
    if fileName[-13:]=="_tracks_1.csv":
        paramFile = FOLDER+"/"+fileName[:-13]
    elif fileName[-11:]=="_tracks.csv":
        paramFile = FOLDER+"/"+fileName[:-11]


    if paramFile != None:
        OUTPUT_FILE = paramFile+"_mixing.csv"
        # os.remove(OUTPUT_FILE)
        # Calculate mixing index for every time step and realisations and save it in OUTPUT_FILE
        write_data_for_all_realisations("mixingIdx", paramFile, OUTPUT_FILE, NEIGHBOUR_DISTANCE)



# #Parameters
# neighbourDistance = 1.05

# for PATH in PATHS:
#     green, red, params = load_simulation_data(PATH)

#     # Find last time step
#     assert(green.find_t_max()==red.find_t_max()), "Tracks don't finish at the same time"
#     t = green.find_t_max()

#     #Get mixing index
#     mixingIdx, mixing = calculate_mixing_index([green,red], t, neighbourDistance = neighbourDistance)

#     print(f"Mixing index: {mixingIdx}")

#     #Write mixing index into file
#     paramsToWrite = {"Pe":params["Pe"], "areaFraction":params["areaFraction"]}
#     strIdx = PATH.rfind("/")
#     FILE_PATH = PATH[:strIdx+1]+"mixingIndex"
#     write_mixing_index(paramsToWrite, mixingIdx, FILE_PATH)

