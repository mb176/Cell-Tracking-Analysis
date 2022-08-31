# Overview

This code is used to analyse and visualise cell track data. It can process the .csv files produced by agent_simulation as well as the experimental data from .xml files. It can compute a variety of observables from the tracks, namely MSD, toruosity, radial distribution functions and mixing index. The functions computing these features are all found in analysis_library.py, and the rest of the files are there to set up to use these functions to calculate and plot specific observables. Files whose names end with _tracks.py analyse experimental tracks, files with the _simulation.py ending apply to simulation results. 

# File list

- **analysis_library.py**: Defines the experiment class, which serves to hold and analyse simulation or experimental data. It can read tracks from .csv and .xml files as well as parameter files. The experiment class cannot distinguish between multiple types of particles, so when we are looking for green-red interactions, we have to load two experiments and use some of the additional functions outside the class. In addition the file contains some stand alone functions that can be used for visualisation. 
- **animation_simulation.py**: Creates an movie showing the simulation tracks.
- **cluster_analysis_simulation.py**:  Calculates the cluster sizes and plots a histogram.
- **cluster_analysis_tracks.py**: Same as before but for the experimental data.
- **final_snapshot_simulation**: Plots and saves a picture of the final arrangement of particles.
- **mixing_index_simulation.py**: Takes a paramter file as input, which is given either by modifying the file or as a command-line argument. It finds the tracks belonging to the parameter file (assuming they are in the same folder) and splits them into green and red particles, which are stored in two different experiments. Then it calculates their mixing index and writes the result with some parameters into a file named "\<parameter_file>_mixing_index.
- **mixing_index_tracks.py**: Calculates the mixing index for the experimental data and plots it.
- **MSD_simulation.py**: Calculates and plots the mean squared displacement for red and green particles and compares them to theoretical predictions.
- **plot_mixing_index_simulation.py**: This is to plot the mixing index at different points in parameter space. To this end it reads the file created by previous executions of mixing_index_simulation.py and plots the data points from there.
- **RDF_simulation.py**: Loads red and green track from the given simulation. Calculates tortuosity (with or without Dun method) as well as the radial distribution function and plots those.
- **test_pytest.py**: Holds unit tests of analysis_library.py. Not fully developed.
- **track_analysis.py**: Similar to RDF_simulation.py but for experimental data: Computes tortuosity and radial density function from tracks and plots the results.
