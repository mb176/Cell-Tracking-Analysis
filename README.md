# Overview

This code is used to analyse and visualise cell track data. It can process the .csv files produced by agent_simulation as well as the experimental data from .xml files. It can compute a variety of observables from the tracks, namely MSD, toruosity, radial distribution functions and mixing index. The functions computing these features are all found in analysis_library.py, and the rest of the files are there to set up to use these functions to calculate and plot specific observables. Some observables, namely mixing index and clustering, require you to analyse the tracks with a "write"-file first, which prints the observable into a file, and then a "plot"-file reads the results and turns them into a plot.
Files whose names end with _exp.py analyse experimental tracks, files with the _sim.py ending apply to simulation results. 

# File list

- **analysis_library.py**: Defines the experiment class, which serves to hold and analyse simulation or experimental data. It can read tracks from .csv and .xml files as well as parameter files. The experiment class cannot distinguish between multiple types of particles, so when we are looking for green-red interactions, we have to load two experiments and use some of the additional functions outside the class. In addition the file contains some stand alone functions that can be used for visualisation. 
- **animation_sim.py**: Creates an movie showing particles at each time step for a given simulation.
- **cluster_analysis_exp.py**: Calculates the cluster sizes and plots a histogram for a given experiment. 
- **clustering_validation.py**: This file checks wether the cluster-finding algorithms work as intended by producing 
a plot where the particles are colored based on cluster.
- **development**: This file exists to test and play around with the functions in analysis_library.py.
- **mixing_index_exp.py**: Calculates the mixing index for the experimental data and plots it.
- **MSD_sim.py**: Calculates and plots the mean squared displacement for red and green particles and compares them to theoretical predictions.
- **plot_clustering_histogram_sim.py**: Reads the clustering.csv file for a given simulation and plots the cluster size histogram for red and green particles. Also checks for statistically relevant differences between the distributions.
- **plot_clustering_over_time_sim.py**: Reads all clustering.csv files in the folder and plots the average maximum cluster size over time for all simulations in one graph (average is over multiple realisations of the same simulation). 
- **plot_last_frame_sim.py**: Plots and saves a picture of the final arrangement of particles.
- **plot_mixing_index_sim.py**: Reads all clustering.csv files in the folder and plots the average mixing index over time for all simulations in one graph (average is over multiple realisations of the same simulation). 
- **RDF_sim.py**: Loads red and green track from the given simulation. Calculates tortuosity (with or without Dun method) as well as the radial distribution function and plots those.
- **track_analysis_exp.py**: Similar to RDF_simulation.py but for experimental data: Computes tortuosity and radial density function from tracks and plots the results.
- **write_clustering_sim.py**: Finds the clustering for all realisations of a given simulation and writes it into a file named "*_clustering.csv".
- **write_mixing_index_sim.py**: Finds the mixing index for all realisations of a given simulation and writes it into a file named "*_mixing.csv".
