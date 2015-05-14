# longscore: longitudinal synapse dynamics data visualization. works with binary presence matrices and supports output from h_imstack2/a loading


# dependencies: glob, numpy, matplotlib.pyplot, seaborn, pandas, scipy.io,scipy.stats
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.io as sio
import scipy.stats as stat
import math as math
import matplotlib as mpl
# DendriticBranch is object for each fov with full timecourse data for spine and shaft densities and dynamics
# bd = by day

print('longscore module loading...')
print('longscore module loaded.')

# Plot and Calc Ranksum p-Value Betwixt Two Vectors


class Dendrite(object):
    def __init__(self,fname):
        # Extract Information from .mat File
        aa = sio.loadmat(fname)
        sad = aa['spineAnalysisData']
        self.total_presence = sad['presenceMatrix'][0,0].astype(int)
        is_spine = sad['isSpineMatrix'][0,0]
        self.dendrite_id = fname
        self.dendrite_length = float(sad['segment_length'][0,0][0,0])

        # Generate Spine and Shaft Presence Matrices
        find_spines = np.argwhere(is_spine)
        self.spine_index = np.unique(find_spines[:,0])
        self.spine_presence = self.total_presence[self.spine_index]
        self.shaft_presence = np.delete(self.total_presence,self.spine_index,0)

        # Calculate Total, Shaft, Spine Numbers and Densities
        self.num_total_bd = self.total_presence.sum(0)
        self.num_spine_bd = self.spine_presence.sum(0)
        self.num_shaft_bd = self.shaft_presence.sum(0)
        self.total_density_bd = self.num_total_bd/self.dendrite_length
        self.spine_density_bd = self.num_spine_bd/self.dendrite_length
        self.shaft_density_bd = self.num_shaft_bd/self.dendrite_length

        ### Calculate Number of Total, Shaft, Spine Dynamics By Day (Additions, Eliminations, Additions + Eliminations)
        # Total
        (self.num_total_dynamic_events_bd,self.num_total_additions_bd,
            self.num_total_eliminations_bd) = calculateSynapseDynamics(self.total_presence)

        # Spine
        (self.num_spine_dynamic_events_bd,self.num_spine_additions_bd,
            self.num_spine_eliminations_bd) = calculateSynapseDynamics(self.spine_presence)

        # Shaft
        (self.num_shaft_dynamic_events_bd,self.num_shaft_additions_bd,
            self.num_shaft_eliminations_bd) = calculateSynapseDynamics(self.shaft_presence)


        #  Calculate Densities of Total, Shaft, and Spine Dynamic Events (Additions, Eliminations, Additions+Eliminations)
        self.total_dynamic_events_density_bd = self.num_total_dynamic_events_bd/self.dendrite_length
        self.total_additions_density_bd = self.num_total_additions_bd/self.dendrite_length
        self.total_eliminations_density_bd = self.num_total_eliminations_bd/self.dendrite_length

        self.spine_dynamic_events_density_bd = self.num_spine_dynamic_events_bd/self.dendrite_length
        self.spine_additions_density_bd = self.num_spine_additions_bd/self.dendrite_length
        self.spine_eliminations_density_bd = self.num_spine_eliminations_bd/self.dendrite_length

        self.shaft_dynamic_events_density_bd = self.num_shaft_dynamic_events_bd/self.dendrite_length
        self.shaft_additions_density_bd = self.num_shaft_additions_bd/self.dendrite_length
        self.shaft_eliminations_density_bd = self.num_shaft_eliminations_bd/self.dendrite_length





class Genotype(object):
    def __init__(self,genomatfiles):
        ### Init Genotype Variables
        # Mean Synapse Density
        self.mean_total_density = []
        self.mean_spine_density = []
        self.mean_shaft_density = []

        # Mean Num Total Dynamics
        self.mean_num_total_eliminations = []
        self.mean_num_total_additions = []
        self.mean_num_total_dynamic_events = []

        # Mean Num Spine Dynamics
        self.mean_num_spine_eliminations = []
        self.mean_num_spine_additions = []
        self.mean_num_spine_dynamic_events = []

        # Mean Num Shaft Dynamics
        self.mean_num_shaft_eliminations = []
        self.mean_num_shaft_additions = []
        self.mean_num_shaft_dynamic_events = []

        # Mean Density Total Dynamics
        self.mean_total_dynamic_events_density = []
        self.mean_total_additions_density = []
        self.mean_total_eliminations_density = []

        # Mean Density Spine Dynamics
        self.mean_spine_additions_density = []
        self.mean_spine_eliminations_density = []
        self.mean_spine_dynamic_events_density = []

        # Mean Density Shaft Dynamics
        self.mean_shaft_dynamic_events_density = []
        self.mean_shaft_additions_denisty = []
        self.mean_shaft_eliminations_density = []



        for i in range(len(genomatfiles)):
            aa = Dendrite(genomatfiles[i])

            # Same Order as Above
            self.mean_total_density.append(np.mean(aa.total_density_bd))
            self.mean_spine_density.append(np.mean(aa.spine_density_bd))
            self.mean_shaft_density.append(np.mean(aa.shaft_density_bd))

            self.mean_num_total_eliminations.append(np.mean(aa.num_total_eliminations_bd))
            self.mean_num_total_additions.append(np.mean(aa.num_total_additions_bd))
            self.mean_num_total_dynamic_events.append(np.mean(aa.num_total_dynamic_events_bd))

            self.mean_num_spine_eliminations.append(np.mean(aa.num_spine_eliminations_bd))
            self.mean_num_spine_additions.append(np.mean(aa.num_spine_additions_bd))
            self.mean_num_spine_dynamic_events.append(np.mean(aa.num_spine_dynamic_events_bd))

            self.mean_num_shaft_eliminations.append(np.mean(aa.num_shaft_eliminations_bd))
            self.mean_num_shaft_additions.append(np.mean(aa.num_shaft_additions_bd))
            self.mean_num_shaft_dynamic_events.append(np.mean(aa.num_shaft_dynamic_events_bd))

            self.mean_total_dynamic_events_density.append(np.mean(aa.total_dynamic_events_density_bd))
            self.mean_spine_dynamic_events_density.append(np.mean(aa.spine_dynamic_events_density_bd))
            self.mean_shaft_dynamic_events_density.append(np.mean(aa.shaft_dynamic_events_density_bd))


# To be called once inside directory with all h_imstack files
class Experiment(object):
    def __init__(self):

        plt.ion()

        # Dictionary Containing Key = Animal Number and Value = String Containing Genotype
        self.geno_dict = {'2193': 'PV','2196':'PV','2374':'VIP'}

        # Load All Matfiles
        self.matfiles = np.array(glob.glob('*.mat'))

        # Init genomatfiles
        self.vipmatfiles = []
        self.pvmatfiles = []
        self.iumatfiles = []


        # Create vip-,pv-,iu-matfiles lists of filenames
        for i in range(len(self.matfiles)):
            if self.geno_dict[self.matfiles[i][:4]] == 'VIP':
                self.vipmatfiles.append(self.matfiles[i])
            elif self.geno_dict[self.matfiles[i][:4]] == 'PV':
                self.pvmatfiles.append(self.matfiles[i])
            elif self.geno_dict[self.matfiles[i][:4]] == 'IU':
                self.iumatfiles.append(self.matfiles[i])

        # Load Genotypes
        self.vip = Genotype(self.vipmatfiles)
        self.pv = Genotype(self.pvmatfiles)
        self.iu = Genotype(self.iumatfiles)

## MODULE METHODS ##

# Calculates Ranksum p-Value and Barplots Two Vectors
def processTwoVectors(vect_one_title,vect_two_title,vect_one,vect_two):
    mean_vect_one = np.mean(vect_one)
    std_vect_one = np.std(vect_one)/math.sqrt(len(vect_one))
    mean_vect_two = np.mean(vect_two)
    std_vect_two = np.std(vect_two)/math.sqrt(len(vect_two))
    bb = stat.ranksums(vect_one,vect_two)
    print(bb)


    N = 2
    means = (mean_vect_one,mean_vect_two)
    stds = (std_vect_one,std_vect_two)
    ind = np.arange(N)
    width = 0.5

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind,means,width,color='c',yerr=stds)
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels((vect_one_title,vect_two_title))

    font = {'family': 'normal', 'weight': 'bold', 'size': 22}
    mpl.rc('font',**font)
    plt.show()

#! Takes Binary Matrix, Returns (Dynamic,Additions,Eliminations)
def calculateSynapseDynamics(binary_presence_matrix):
    diff_matrix = np.diff(binary_presence_matrix)
    dynamic_matrix = np.absolute(diff_matrix)
    additions_matrix = (diff_matrix>0).astype(int)
    eliminations_matrix = (diff_matrix<0).astype(int)
    dynamic_events_bd = dynamic_matrix.sum(0)
    additions_bd = additions_matrix.sum(0)
    eliminations_bd = eliminations_matrix.sum(0)
    return (dynamic_events_bd,additions_bd,eliminations_bd)
