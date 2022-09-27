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
    PATH = sys.argv[1]
else:
    PATH = "/home/marius/PhD/CellMotility/agent_simulation/output/variedPersistence/mixingIndex"
    
xValues, yValues, mixingIndices, names = read_mixing_index(PATH)
fig, ax = plt.subplots(1,1)
im = ax.scatter(xValues, yValues, s=200, c=mixingIndices, cmap='Blues')
ax.set_xlabel(names[0])
ax.set_ylabel(names[1])
fig.colorbar(im)
plt.savefig(PATH+"_phasediagram.png")