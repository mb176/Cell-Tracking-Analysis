import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # This avoids crashes on the math cluster from using to many resources
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
from scipy import stats
from matplotlib import animation
import copy
import sys
from os.path import exists
from analysis_library import *


if(len(sys.argv)==2):
    #Parameter file given externally:
    parameterFile = sys.argv[1]
else:
    #Give parameter file manually
    parameterFile = "/home/marius/PhD/CellMotility/agent_simulation/output_delayed_CIL/test/test"

""" This file goes through the realisations belonging to the parameterFile and accumulates a high 
dimensional histogram that includes self-propulsion angles. The histogram is then written to 
paramterFile_histogram.csv. It doesn't consider phi1, phi2 and phiR independently, but only the
combinations theta1 = phiR-phi1  and theta2 = phi2-phi1
"""


# Defining the histogram
rBinWidth = 0.5          # Bin size for the radial distance
theta1BinWidth = np.pi/12  # Bin size for the angle between the vector connecting two particles and the x axis
theta2BinWidth = np.pi/12  # Bin size for the angle of the self-propulsion velocity of particle 1
times = [0,11,14,17,20]  # All the time slices that are included in the histogram
maxRealisations = 100   # How many are at most included in the histogram

sim = experiment("sim")
sim.periodic = True
params = sim.read_parameter_file(parameterFile)
maxLength = sim.set_max_distance()

nDataPoints = 0
nRBins = int(np.ceil(maxLength/rBinWidth)) #Only positive values
nTheta1Bins = int(np.ceil(2*np.pi/theta1BinWidth)) 
nTheta2Bins = int(np.ceil(2*np.pi/theta2BinWidth)) 

histogram = np.zeros((nRBins, nTheta1Bins, nTheta2Bins, 3 , len(times))) 
# Distance, theta1 = phiR-phi1, theta2 = phi2-phi1, typeIdx, time
# type has: 0 (both red), 1 (both green), 2 (green-red pair)


n=1
while(exists(parameterFile+f"_tracks_{n}.csv") and n <= maxRealisations):
    print(f"Scanning realisation number {n}\n")
    measurementTimes, colorMeasurements, positionMeasurements = readTracksFile(parameterFile+f"_tracks_{n}.csv")
    velocityMeasurements = sim.read_velocity_csv(parameterFile+f"_velocities_tracks_{n}.csv")
    n+=1
    
    for timeIdx, time in enumerate(measurementTimes):
        if time in times:
            print(f"Scanning time slice t={time}\n")
            for pIdx1 in range(int(params["nParticles"])):
                for pIdx2 in range(pIdx1):
                    # Distance & angle
                    r, phiR = sim.distance(np.array(positionMeasurements[timeIdx][pIdx2])-np.array(positionMeasurements[timeIdx][pIdx1]),angle=True)
                    # Velocity angles
                    phi1 = velocityMeasurements[pIdx1][timeIdx]
                    phi2 = velocityMeasurements[pIdx2][timeIdx]
                    if phi1<0: phi1+= 2*np.pi
                    if phi2<0: phi2+= 2*np.pi
                    theta1 = phiR - phi1
                    theta2 = phi2 - phi1
                    if theta1<0: theta1+= 2*np.pi
                    if theta2<0: theta2+= 2*np.pi
                    # Type of pair
                    if(colorMeasurements[timeIdx][pIdx1]==" red" and colorMeasurements[timeIdx][pIdx2]==" red"): 
                        typeIdx = 0
                    elif(colorMeasurements[timeIdx][pIdx1]!=" red" and colorMeasurements[timeIdx][pIdx2]!=" red"):
                        typeIdx = 1
                    else:
                        typeIdx = 2

                    # Type: 0 - red&red 1 - green&green 2 - red&green
                    # Time 
                    tIdx = times.index(time)

                    # Add the connection from each particle's perspective
                    try: # Very rarely you can get out of bound indices thorugh rounding error
                        histogram[int(r/rBinWidth), 
                                int(theta1/theta1BinWidth), 
                                int(theta2/theta2BinWidth),
                                typeIdx, 
                                tIdx] += 1
                    except:
                        continue

    
cpy = histogram 
# Save the histogram to file
outputFile = parameterFile+"_histogram_delta.csv"
with open(outputFile, 'w') as f:
    f.write("Frequency, distance, theta1 = phiR-phi1, theta2 = phi2-phi1, type of pair, time slice \n")
    f.write(f"binWidths: {rBinWidth} {theta1BinWidth} {theta2BinWidth} {len(times)} \n")
    for rIdx in range(nRBins):
            for theta1Idx in range(nTheta1Bins):
                for theta2Idx in range(nTheta2Bins):
                    for typeIdx in range(3):
                        for timeIdx, time in enumerate(times):
                            f.write(f"{histogram[rIdx, theta1Idx, theta2Idx, typeIdx, timeIdx]} {rIdx} {theta1Idx} {theta2Idx} {typeIdx} {time} \n")
