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
#python celluloid
#You can use np.arrays with named vectors, just like pandas!
# rain_drops = np.zeros(n_drops, dtype=[('position', float, (2,)),
#                                       ('size',     float),
#                                       ('growth',   float),
#                                       ('color',    float, (4,))])

# # Initialize the raindrops in random positions and with
# # random growth rates.
# rain_drops['position'] = np.random.uniform(0, 1, (n_drops, 2))

def readTracksFile(filename):
    """This is similar to read_csv in analysis library, but has a different output structure and is
    only used here"""
    measurementTimes = []
    positionMeasurements = []
    colorMeasurements = []
    with open(filename, mode="r") as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            measurementTimes.append(float(lines[0])) 
            positionMeasurements.append([])
            colorMeasurements.append([])
            nParticles = int((len(lines)-1)/3)
            for particleIdx in range(nParticles):
                colorMeasurements[-1].append(lines[1+3*particleIdx])
                positionMeasurements[-1].append([float(lines[1+3*particleIdx+1]), #x value
                                              float(lines[1+3*particleIdx+2])])   #y value

    return measurementTimes, colorMeasurements, np.array(positionMeasurements)
     
def getColors(colorMeasurement): ##ToDo: Get right colors
    c = []
    for color in colorMeasurement:
        if color==" red":
            c.append('red')
        elif color==" green":
            c.append('green')
        elif color==" greenPlus":
            c.append('purple')
        else:
            assert False, "Invalid color value"
    return c


#SOURCE_FOLDER = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', "output"))+"/" #wtf python
#NAME = "test_parameters_tracks"

#Source file
PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters"]
# PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity/LowDensity",
#         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity/HighDensity",
#         "/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl/HighDensityControl"]

for PATH in PATHS:
    #Get parameters
    sim = experiment("test")
    params = sim.read_parameter_file(PATH)
    xLength =params["Length"] 
    yLength =params["Length"]
    R = 0.5 #Radius of the cells since we set r0=1

    #Get tracks
    measurementTimes, colorMeasurements, positionMeasurements = readTracksFile(PATH+"_tracks.csv")
    positionMeasurements = np.array(positionMeasurements)


    #Set up figure 
    fig,ax = plt.subplots()
    plt.axis([0, xLength, 0, yLength])
    plt.gca().set_aspect('equal', adjustable='box') #Makes both axis have some scaling while keeping the set limits

    #Record the animation
    camera = Camera(fig)
    for frameIdx in range(len(measurementTimes)):
        c = getColors(colorMeasurements[frameIdx])
        for particleIdx in range(len(positionMeasurements[0])):
            circle = plt.Circle((positionMeasurements[frameIdx][particleIdx][0],
                                        positionMeasurements[frameIdx][particleIdx][1]), 
                                        radius=R, linewidth=0, color = c[particleIdx])
            plt.gca().add_artist(circle)
            plt.gca().set_title("Time =%f"%measurementTimes[frameIdx])
        camera.snap()
        #Delete all artists before next frame
        for artist in plt.gca().lines + plt.gca().collections:
            artist.remove()
    animation = camera.animate()
    animation.save(PATH+".mov")

