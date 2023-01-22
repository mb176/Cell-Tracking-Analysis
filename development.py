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

# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/marius/PhD/CellMotility/analysis/')
from analysis_library import *

"""
This file exists to test and play around with the 
functions in analysis_library.py.
"""

def read_histogram_file(parameterFile):
    """Reads a previously written file containing a multidimensional histogram. It assumes that the histogram is 
    in the same folder as the parameterFile and is named parameterFile+"_histogram.csv"""
    sim = experiment("sim")
    sim.periodic = True
    params = sim.read_parameter_file(parameterFile)
    maxLength = sim.set_max_distance()
    with open(parameterFile+"_histogram.csv", 'r') as f:
        for line in f.readlines():
            if(line[0]=="Frequency"): continue
            elif(line[0]=="binWidth"):
                rBinWidth = float(line[1])
                phiRBinWidth = float(line[2])
                phi1BinWidth = float(line[3])
                phi2BinWidth = float(line[4])
                nTimes = float(line[5])
                nRBins = int(np.ceil(maxLength/rBinWidth)) #Only positive values
                nPhiRBins = int(np.ceil(2*np.pi/phiRBinWidth)) 
                nPhi1Bins = int(np.ceil(2*np.pi/phi1BinWidth)) 
                nPhi2Bins = int(np.ceil(2*np.pi/phi2BinWidth))
                histogram = np.zeros((nRBins, nPhiRBins, nPhi1Bins, nPhi2Bins, 3)) 

            else:
                histogram[  int(float(line[1])/rBinWidth), 
                            int(float(line[2])/phiRBinWidth), 
                            int(float(line[3])/phi1BinWidth),
                            int(float(line[4])/phi2BinWidth),
                            int(line[5])] = int(line[0])


if(len(sys.argv)==2):
    #Parameter file given externally:
    PATHS = [sys.argv[1]]
else:
    #Give parameter file manually
    PATHS = ["/home/marius/PhD/CellMotility/agent_simulation/output_delayed_CIL/test/test"
            ]

parameterFile = PATHS[0]
rBinWidth = 0.5
phiRBinWidth = np.pi/12
phi1BinWidth = np.pi/12
phi2BinWidth = np.pi/12
times = [11, 13, 15, 17, 19]
timeAvr = True
pairType = "all" #"red" "green", "mixed", "all"
"""Reads a previously written file containing a multidimensional histogram. It assumes that the histogram is 
in the same folder as the parameterFile and is named parameterFile+"_histogram.csv"""
sim = experiment("sim")
sim.periodic = True
params = sim.read_parameter_file(parameterFile)
maxLength = sim.set_max_distance()

sim = experiment("sim")
sim.periodic = True
params = sim.read_parameter_file(parameterFile)
maxLength = sim.set_max_distance()

nDataPoints = 0
nRBins = int(np.ceil(maxLength/rBinWidth)) #Only positive values
nPhiRBins = int(np.ceil(2*np.pi/phiRBinWidth)) 
nPhi1Bins = int(np.ceil(2*np.pi/phi1BinWidth)) 
nPhi2Bins = int(np.ceil(2*np.pi/phi2BinWidth))

histogram = np.zeros((nRBins, nPhiRBins, nPhi1Bins, nPhi2Bins, 3 , len(times))) 
# Distance, angle between particles, velcity angle 1, velocity angle 2, type , time slice, 
# type has: 0 (both red), 1 (both green), 2 (green-red pair)


n=1
while(exists(parameterFile+f"_tracks_{n}.csv") and n<6):
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
                    # Type of pair
                    if(colorMeasurements[timeIdx][pIdx1]==" red" and colorMeasurements[timeIdx][pIdx1]==" red"): 
                        typeIdx = 0
                    elif(colorMeasurements[timeIdx][pIdx1]!=" red" and colorMeasurements[timeIdx][pIdx1]!=" red"):
                        typeIdx = 1
                    else:
                        typeIdx = 2
                    # Time 
                    tIdx = times.index(time)

                    # Add the connection from each particle's perspective
                    histogram[int(r/rBinWidth), 
                              int(phiR/phiRBinWidth), 
                              int(phi1/phi1BinWidth),
                              int(phi2/phi2BinWidth),
                              typeIdx, 
                              tIdx] += 1
                    if(phiR >= np.pi): phiRPrime = phiR - np.pi
                    else: phiRPrime = phiR + np.pi
                    histogram[int(r/rBinWidth), 
                              int(phiRPrime/phiRBinWidth), 
                              int(phi2/phi2BinWidth),
                              int(phi1/phi1BinWidth),
                              typeIdx, 
                              tIdx] += 1
                    assert(phiR>=0)
                    assert(phiRPrime >= 0)
                    nDataPoints += 2
    
cpy = histogram 
print(histogram.max())

# Save the histogram to file
outputFile = parameterFile+"_histogram.csv"
with open(outputFile, 'w') as f:
    f.write("Frequency, distance, angle between particles, velocity angle 1, velocity angle 2, type of pair, time slice \n")
    f.write(f"binWidths: {rBinWidth} {phiRBinWidth} {phi1BinWidth} {phi2BinWidth} {len(times)} \n")
    for rIdx in range(nRBins):
        for phiRIdx in range(nPhiRBins):
            for phi1Idx in range(nPhi1Bins):
                for phi2Idx in range(nPhi2Bins):
                    for typeIdx in range(3):
                        for timeIdx, time in enumerate(times):
                            f.write(f"{histogram[rIdx, phiRIdx, phi1Idx, phi2Idx, typeIdx, timeIdx]} {rIdx} {phiRIdx} {phi1Idx} {phi2Idx} {typeIdx} {timeIdx} {time} \n")
                            print(f"{histogram[rIdx, phiRIdx, phi1Idx, phi2Idx, typeIdx, timeIdx]} {rIdx} {phiRIdx} {phi1Idx} {phi2Idx} {typeIdx} {timeIdx} {time} \n")
# Reference histogram
minimum = 0
nReferencePoints = 5000
referenceHistogram = sim.generate_reference_histogram(nReferencePoints,nRBins, nAngleBins = nPhiRBins)


# Cutoff all values that are in the corners 
cutoff = int(maxLength/np.sqrt(2)/rBinWidth)
histogram = histogram[0:cutoff]
referenceHistogram = referenceHistogram[0:cutoff]

# Marginalize the distribution (needs to go from last to first index)
if(timeAvr==True):
    histogram = np.sum(histogram, 5)
if(pairType=="all"):
    histogram = np.sum(histogram, 4)
if(nPhi2Bins==1):
    histogram = np.sum(histogram, 3)
if(nPhi1Bins==1):
    histogram = np.sum(histogram, 2)
if(nPhiRBins==1):
    histogram = np.sum(histogram, 1)


#Normalise the radial density
nReferenceValues = (nReferencePoints-1)*(nReferencePoints)/2
for idx, value in enumerate(histogram[0,:]):
    nValues = np.sum(histogram[:,idx])
    if(nValues>0):
        histogram[:,idx] = nReferenceValues/nValues*histogram[:,idx]/referenceHistogram

#Normalise the angular density
if nPhiRBins > 1:
    referenceHistogram *=1/nReferenceValues
    for idx, value in enumerate(histogram[:,0]):
        histogram[idx,:] = histogram[idx,:]/referenceAngleHistogram

#Do both at same time:
if nPhiRBins > 1:
    ref = np.outer(referenceHistogram, referenceAngleHistogram)
    histogram = nReferenceValues/nValues*histogram /ref

# Create bin
rBins = np.linspace(0, maxLength, cutoff+1)[:-1] + rBinWidth #Center of each bin
phi1Bins = np.linspace(0, 2*np.pi, nPhi1Bins+1)[:-1] + phi1BinWidth

# 1D plot
fig, ax = plt.subplots(1,1) 
ax.plot(rBins, histogram)

# 2D plot
fig3D = plt.figure()
ax3D = plt.axes(projection="3d")
X = np.outer(rBins, np.ones(1))
Y = np.outer(np.ones(1), phi1Bins)
ax3D.plot_surface(X, Y, histogram)

