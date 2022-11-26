import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
from analysis_library import *
plt.style.use('seaborn-whitegrid')
plt.close("all")

"""
This file checks wether the cluster-finding algorithms work as intended by producing 
a plot where the particles are colored based on cluster.
"""

def plot_cluster(axes, experiments, time, radius, neighbourDistance):
        global colorIdx 

        # Get points
        points=[]
        color=[]
        for expIdx, exp in enumerate(experiments):
                for pIdx, track in enumerate(exp.tracks):
                        tIdx, = np.where(np.array(track)[:,0]==time)

                        if len(tIdx)>0: #Does the tracked particle exist at t?
                                tIdx = tIdx[0] #first occurence
                                color.append(expIdx)
                                points.append(np.array(track[tIdx][1:]))


        #Find cluster
        if neighbourDistance==None:
                G = get_voronoi_based_neighbourhood_graph(points, color, experiments[0])
        else:
                G = get_distance_based_neighbourhood_graph(points, color, experiments[0], neighbourDistance)
        

        xLength = green.x_max
        yLength = green.y_max

        # Plot cluster
        for comp in nx.connected_components(G):
                for idx in comp: 
                        assert(color[idx]==color[list(comp)[0]]), "Clusters have to be unicolored!"
                        circle = plt.Circle((points[idx][0],
                                        points[idx][1]), 
                                        radius=radius, linewidth = 0.1 ,alpha = TRANSPARENCY)
                        circle.set_facecolor(f"C{colorIdx}")
                        if color[idx]==0:
                                circle.set_edgecolor("red")   
                        elif color[idx]==1:
                                circle.set_edgecolor("green")
                        else:
                                assert(False, "Invalid Color")
                                       
                
                        axes.add_artist(circle)
                colorIdx += 1 
        
        axes.axis([0, xLength, 0, yLength])
        axes.set_aspect('equal', adjustable='box') #Makes both axis have some scaling while keeping the set limits

#Give parameter file manually
PATH = "/home/marius/PhD/CellMotility/agent_simulation/new_output/turnAround_persistence/areaFraction_0.7_Pe_80"
       


#Parameters
NEIGHBOUR_DISTANCE = 1.05
TRANSPARENCY = 0.8
R = 0.5


green, red, params = load_simulation_data(PATH, trackIdx="_1")

# Find last time step
assert(green.find_t_max()==red.find_t_max()), "Tracks don't finish at the same time"
t = 0 #green.find_t_max()


# Red and green separate
fig2, axes2 = plt.subplots(1,2)  
colorIdx = 0
plot_cluster(axes2[0], [green,red], t, R, NEIGHBOUR_DISTANCE)
axes2[0].set_title(f"Distance (d={NEIGHBOUR_DISTANCE})")

colorIdx = 0
plot_cluster(axes2[1], [green, red], t, R, None)
axes2[1].set_title("Voronoi")

fig2.savefig(PATH+'_initial_cluster.png', format='png',dpi=500)

        # nx.draw_networkx(G)



plt.close('all')