import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
import time
from numba import jit
from scipy import spatial, ndimage

@jit(nopython=True)
def generate_uniform_histogram(n_particles, n_bins, x_max,y_max, seed=None):
        '''Generates a radial density histogram for uniform distributed particles in a box of 
        size x_max y_max. This serves as a normalisation for the measured density.
        :n_particles Number of particles that are generated
        :n_bins Number of bins to be used. The bins are evenly distributed between 0 and max_distance
            where max_distance is calculated from x_max and y_max
        :x_max, y_max Dimensions of the rectangle for which we generate the density
        :rng Random number generator to be used
        '''

        #Numba compatible random number genration
        if seed!=None:
            np.random.seed(seed)

        #Generate all points
        max_distance = np.sqrt(x_max**2+y_max**2)
        x_values = np.random.random((n_particles))*x_max
        y_values = np.random.random((n_particles))*y_max
        points = np.stack((x_values,y_values),axis=1)


        # distances = spatial.distance.pdist(points)
        # reference_histogram = ndimage.histogram(distances, 0, max_distance,n_bins)

        #Prepare histogram
        bins = np.linspace(0, max_distance, n_bins+1)[:-1] #left edge of every bin
        histogram = np.zeros(n_bins)

        #calculate distances between all points
        for n in range(n_particles):
            for m in range(n):
                distance = np.linalg.norm(points[n]-points[m])
                bin_idx = int(np.floor(distance/max_distance*n_bins))
                histogram[bin_idx]+=1
        
        reference_histogram = histogram

        return reference_histogram

class experiment:
    """Holds tracks of the experiment, as well as any calculated observables.
    :name Name of the experiment
    :tracks Array of numpy arrays containing the times and coordinates of each particle
    :tmax number of timesteps
    :x_max Biggest x coordinate found in a track
    :y_max Biggest y coordinate found in a track 
    :distance_travlled Array with one entry for every track, contains the total distance a 
        particle has travelled on the track
    :tortuosity Array with one entry for every track, contains the total normalised
        difference between distance_travelled and the total displacement
    :rng Random number generator used to create histogram normalisaiton by generating particles with 
        uniform density
    :uniform_histogram Stores histogram of radial density of uniform particles. Used
        for normalisation in plot_radial_density
    """
    def __init__(self,name,seed=None):
        self.name = name
        self.tmax = 0
        self.tracks = [] 
        self.x_max = 0 
        self.y_max = 0
        self.reference_histogram = np.array([None])
        self.seed=seed
        if seed==None:
            self.rng = np.random.default_rng()
        else:
            self.rng = np.random.default_rng(self.seed)
    
    def read_xml(self,location, min_length):
        """Reads the xml file in location and saves the tracks in it into location
        :location Path to the file
        :min_length All tracks shorter than that are excluded
        """
        tic = time.perf_counter()

        #Read xml file into beautiful soup:
        with open(location, 'r') as f:
             data = f.read()
        Bs_data = BeautifulSoup(data, "xml")

        #keep track of the biggest x and y values we find

        #Iterate through all tracks in the file
        for track in Bs_data.find_all('particle'):
            #create numpy arrays to contain the track
            nSpots = int(track['nSpots'])
            if nSpots >=min_length:#possible to filter short tracks
                self.tracks.append(np.zeros(shape=(nSpots,3)))
                n=0
                #Iterate through each spot detected for the particle
                for spot in track.find_all('detection'):
                    self.tracks[-1][n][0]=int(spot['t'])
                    self.tracks[-1][n][1]=float(spot['x'])
                    self.tracks[-1][n][2]=float(spot['y'])
                    self.x_max = max(self.x_max,self.tracks[-1][n][1])
                    self.y_max = max(self.y_max,self.tracks[-1][n][2])
                    n+=1
        toc = time.perf_counter()
        print(f"Read xml file in {toc - tic:f} seconds")
        
    
    # def TAMSD(self,track,delta_t):
    #     """Calculates the time averaged mean squared displacement of a single track
    #     with timestep delta_t (given as a time index, not in minutes). 
    #     """
    #     track_duration = len(track)
    #     #Check that track is long enough
    #     assert track_duration > delta_t+1, 'Track duration is shorter than the time interval'
    #     total = 0
    #     for idx in range(track_duration-delta_t):
    #         distance = track[idx][1:]-track[idx+delta_t][1:]
    #         total += np.dot(distance,distance)
    #     return total/(track_duration-delta_t)

    def TAMSD(self,delta_t,min_length):
        """Returns the average TAMSD as well as all individual TAMSDs for all tracks in the experiment
         whose length exceed min_length
        :delta_t time interval (in timesteps, not minutes) for which the TAMASD will be computed
        :min_length minimum number of timesteps a track must have to be used for the average
        """
        assert delta_t < min_length, "The minimum track length must exceed the time average interval"
        nTracks = 0 # counts how many tracks are included/ have the minimum length
        average = 0
        individual_TAMSD = [] 
        for track in self.tracks:
            track_duration = len(track)
            total = 0
            if track_duration >= min_length:
                nTracks += 1
                #Calculate TAMSD for the track
                for idx in range(track_duration-delta_t):
                    distance = track[idx][1:]-track[idx+delta_t][1:]
                    total += np.dot(distance,distance)
                TAMSD= total/(track_duration-delta_t)
                individual_TAMSD.append(TAMSD)
                average += TAMSD
        assert nTracks > 0, 'None of the tracks have the minimum length'
        average = average/nTracks

        return average, individual_TAMSD

    def plot_TAMSD(self, tmax, min_length,axes):
        """Calculates the TAMSD as a function of the time intervall delta_t for all tracks 
        in the experiment that are longer than min_length.
        :tmax sets the biggest delta_t that is considered (in timesteps, not minutes). Has to 
            be an integer
        :min_length is used to filter all tracks that are shorter 
        :axes is the pyplot object into which the graph is drawn
        """
        delta_t_vec = np.linspace(1, tmax, tmax, dtype = int)
        avr_TAMSD = np.zeros(tmax)
        individual_TAMSD = []
        for delta_t in delta_t_vec:
            avr_TAMSD[delta_t-1], TAMSDs = self.TAMSD(delta_t,min_length)
            individual_TAMSD.append(TAMSDs)

        #Plot individual tracks
        individual_TAMSD = np.transpose(np.array(individual_TAMSD))
        for vec in individual_TAMSD:
            axes.loglog(delta_t_vec,vec,linewidth=0.1)    
        axes.set_xlabel(r"$\Delta$ t")
        axes.set_ylabel("TAMSD")
        axes.loglog(delta_t_vec,avr_TAMSD,label='Average', linewidth=2)

        #Plot dashed lines for comparison
        
        #Ballistic motion, starting on the left
        x0 = avr_TAMSD[0]
        x1 = avr_TAMSD[0]*(tmax)**2 
        axes.loglog([delta_t_vec[0],delta_t_vec[-1]],[x0,x1],color='black', 
                    linestyle='dashed', label='Ballistic motion')

        #Normal diffusive motion, 
        y0 = avr_TAMSD[0]
        y1 = avr_TAMSD[0]*(tmax) 
        axes.loglog([delta_t_vec[0],delta_t_vec[-1]],[y0,y1],color='blue', 
                    linestyle='dashed', label='Diffusive motion')

    def calculate_tortuosity_dun_method(self, delta_t, min_length):
        """calculates distance travelled and 
        tortuosity = (distance_travelled-displacement)/distance_travlled for every track and
        saves the values in two arrays, whose indices correspond to those of the tracks
        :delta_t Size of time interval for averaging (Dun method)
        :min_length Minimum size the track must have
        """
        assert delta_t < min_length, "Can't average unless delta_t<min_length"
        nTracks = len(self.tracks)
        self.distance_travelled = np.zeros(nTracks)
        self.tortuosity = np.zeros(nTracks)
        n=0 #track index
        for track in self.tracks:
            track_duration = len(track)
            if track_duration >= min_length:
                tortousity = 0
                n_intervals = int(track_duration//delta_t)
                for idx in range(n_intervals): 
                    #Calculate Tortuosity in every interval
                    straight_line_distance = np.linalg.norm(track[idx*delta_t][1:]-track[(idx+1)*delta_t-1][1:])
                    displacement = 0
                    for spot_idx in range(idx*delta_t,(idx+1)*delta_t-1):
                        displacement += np.linalg.norm(track[spot_idx+1][1:]-track[spot_idx][1:])
                    tortousity += (displacement-straight_line_distance)/displacement
                tortuosity = tortousity/n_intervals

                self.distance_travelled[n]=displacement
                self.tortuosity[n] = tortuosity
                n+=1
        #Extract non-zero entries
        self.distance_travelled = self.distance_travelled[self.tortuosity!=0]
        self.tortuosity = self.tortuosity[self.tortuosity!=0]

    def calculate_tortuosity(self, min_length): #Old version without dun method
        """calculates distance travelled and 
        tortuosity = (distance_travelled-displacement)/distance_travlled for every track and
        saves the values in two arrays, whose indices correspond to those of the tracks
        """
        nTracks = len(self.tracks)
        self.distance_travelled = np.zeros(nTracks)
        self.tortuosity = np.zeros(nTracks)
        n=0 #track index
        for track in self.tracks:
            if len(track)>= min_length:
                displacment = np.linalg.norm(track[0][1:]-track[-1][1:])
                distance_travelled = 0
                for idx in range(len(track)-1):
                    distance_travelled += np.linalg.norm(track[idx+1][1:]-track[idx][1:])
                tortuosity = (distance_travelled-displacment)/distance_travelled

                self.distance_travelled[n]=distance_travelled
                self.tortuosity[n] = tortuosity
                n+=1
        #Extract non-zero entries
        self.distance_travelled = self.distance_travelled[self.tortuosity!=0]
        self.tortuosity = self.tortuosity[self.tortuosity!=0]
        
    def plot_tortuosity(self,axes, min_length, label, color, n_bins=None, method='Normal', delta_t=None):
        """Plots tortuosity either by considering the entire track at once, or by splitting 
        the track into intervals of length delta_t and averaging over those (Dun method).
        """
        if method=='Normal':
            self.calculate_tortuosity(min_length)
        elif method == 'Dun':
            self.calculate_tortuosity_dun_method(delta_t, min_length)
        else:
            raise Exception("Unknown method for computing tortuosity")
        #axes.scatter(self.distance_travelled,self.tortuosity, s=0.4)
        if n_bins == None:
            axes.hist(self.tortuosity, density=True, alpha=0.5, label=label, color=color)
        else: 
            axes.hist(self.tortuosity, density=True, alpha=0.5, label=label, color=color, bins=n_bins)
        axes.set_xlabel(r'Tortuosity = (L-I)/L')
        axes.set_ylabel(r'Probability Density')
        axes.legend()

    def generate_uniform_histogram(self, n_particles, n_bins):
        '''Generates a radial density histogram for uniform distributed particles in a box of 
        size self.x_max self.y_max. This serves as a normalisation for the measured density.
        :n_particles Number of particles that are generated
        :n_bins Number of bins to be used. The bins are evenly distributed between 0 and max_distance
            where max_distance is calculated from self.x_max and self.y_max
        '''
        #Check if we already created the distogram with this number of particles

        # #Generate all points
        # max_distance = np.sqrt(self.x_max**2+self.y_max**2)
        # x_values = self.rng.random((n_particles))*self.x_max
        # y_values = self.rng.random((n_particles))*self.y_max
        # points = np.transpose(np.array([x_values,y_values]))


        # distances = spatial.distance.pdist(points)
        # self.reference_histogram = ndimage.histogram(distances, 0, max_distance,n_bins)

        if (all(self.reference_histogram==None )):  #Check if we already created a reference histogram of this size
            #or np.sum(self.reference_histogram)!=(n_reference_points-1)*(n_reference_points)/2):
            self.reference_histogram = generate_uniform_histogram(n_particles, n_bins, self.x_max,
                                    self.y_max, self.seed)
            

        return self.reference_histogram

    def plot_radial_density(self,axes, t ,n_bins, label, no_plot=False, cutoff_percentange = 10, 
                            n_reference_points = 10000):
        """Goes through all tracks to find spots in the time slice t and finds all distances 
        between them. The distances are assigned to the n_bin histogram bins.
        :axes Histogram will be plotted here
        :t Time at which the correlation function will taken
        :n_bins number of bins in the histogram
        :no_plot Does not plot on axes
        :cutoff_percentage The percentage of bins that will be cutoff at the end.
            The tail end 
        :n_reference_points Number of points used for the reference histogramm. High values 
            tax performance, but smoothen the reference density
        """
        #Prepare histogram
        max_distance = np.sqrt(self.x_max**2+self.y_max**2)
        bins = np.linspace(0, max_distance, n_bins+1)[:-1] #left edge of every bin
        histogram = np.zeros(n_bins)
        n_points = 0

        #Search all tracks for points
        points = []
        for track in self.tracks:
            t_start = track[0][0]
            if t_start <= t and t_start+len(track)>t: #Does the tracked particle exist at t?
                new_point = track[int(t-t_start)][1:]
                n_points += 1
                #Compute all distances to the new point and place them in the histogram
                for old_point in points: 
                    distance = np.linalg.norm(new_point-old_point)
                    bin_idx = int(np.floor(distance/max_distance*n_bins))
                    histogram[bin_idx]+=1
                points.append(new_point)


        #Normalise the density
        n_distances = (n_points-1)*(n_points)/2
        n_reference_distances = (n_reference_points-1)*(n_reference_points)/2
        

        tic = time.perf_counter()
        self.reference_histogram = self.generate_uniform_histogram(n_reference_points,n_bins)
        toc = time.perf_counter()
        print(f"Reference Histogram with {n_reference_points:.0f} particles was computed in {toc - tic:f} seconds")
            
        #Multiply by the ratio of the numbers of distances between points in each histogram
        reference_density = (n_distances/n_reference_distances)*self.reference_histogram
        cutoff = int(np.floor(n_bins*cutoff_percentange/100)) #Cutoff tail end of the histogarm to avoid dividing by 0
        radial_density = np.divide(histogram[:-cutoff],reference_density[:-cutoff])

        #Plot Histogram
        if no_plot==False:
            bin_midpoints = bins+(bins[1]-bins[0])/2
            axes.plot(bin_midpoints[:-cutoff], radial_density, label=label)
            axes.set_xlabel(r'Distance')
            axes.set_ylabel(r'Radial distribution function')
            axes.legend()
        
        return histogram, bins

    

def plot_mixed_particle_radial_density(axes, experiments, t ,n_bins, n_reference_points, no_plot=False, cutoff_percentage=0):
    """Assumes input of 2 experiments
    """
    max_distance = 0
    for exp in experiments:
        max_distance = max(max_distance,np.sqrt(exp.x_max**2+exp.y_max**2))
    bins = np.linspace(0, max_distance, n_bins+1)[:-1] #left edge of every bin
    histogram = np.zeros(n_bins)
    n_points = 0
    
    tic = time.perf_counter()
    position_array = []
    for exp in experiments:
        points = []
        for track in exp.tracks:
            t_start = track[0][0]
            if t_start <= t and t_start+len(track)>t: #Does the tracked particle exist at t?
                points.append(track[int(t-t_start)][1:])
        position_array.append(points) #has spots from both experiments
    points1 = np.array(position_array[0])
    points2 = np.array(position_array[1])
    toc = time.perf_counter()
    print(f"Read out all spots of both experiments in one time slice in {toc - tic:f} seconds")

    # tic = time.perf_counter()
    # for point1 in position_array[0]:
    #     for point2 in position_array[1]:
    #         distance = np.linalg.norm(point1-point2)
    #         bin_idx = int(np.floor(distance/max_distance*n_bins))
    #         histogram[bin_idx]+=1
    # toc = time.perf_counter()
    # print(f"Computed all distances with for loops in {toc - tic:f} seconds")

    tic = time.perf_counter()
    distances = spatial.distance.cdist(points1,points2)
    histogram = ndimage.histogram(distances, 0, max_distance,n_bins)
    toc = time.perf_counter()
    print(f"Computed all distances with scipy in {toc - tic:f} seconds")



    reference_histogram = experiments[0].generate_uniform_histogram(n_reference_points,n_bins)
    #Multiply by the ratio of the numbers of links between points in each histogram
    n_distances = len(position_array[0])*len(position_array[1])
    n_reference_distances = n_reference_points*(n_reference_points-1)
    reference_histogram = (n_distances/n_reference_distances)*reference_histogram
    cutoff = int(np.floor(n_bins*cutoff_percentage/100)) #Cutoff tail end of the histogarm to avoid dividing by 0
    radial_density = np.divide(histogram[:-cutoff],reference_histogram[:-cutoff])

    #Plot Histogram
    if no_plot==False:
        bin_midpoints = bins+(bins[1]/2)
        axes.plot(bin_midpoints[:-cutoff], radial_density, label= 'Red-green RDF')
        axes.set_xlabel(r'Distance')
        axes.set_ylabel(r'Radial distribution function')


