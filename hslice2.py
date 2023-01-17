import math
import time
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

cook_time = 61 # cook time in ms
stack_period = 25 # stack period (duty cycle) in ms
cycles = 38

plot_energy_bounds = [0,10000] # [upper, lower] photon energy axis limits (eV)

energy_bin = 3 # photon energy bins size (eV)
time_bin = 5 # time bins size (ms)


def plot_hslice(ranges):

    data = np.load('/home/pcuser/Desktop/TES-GUI/M3 Lifetime/20221220_0000_hjl.npy')
    
    data = data[(data[:,0] > plot_energy_bounds[0]) & (data[:,0] <= plot_energy_bounds[1])] # truncate data out of energy bounds
    data[:,1] *= 1000
    #data = data[data[:,1] < stack_period*cycles]
    #data[:,1] = data[:,1] % stack_period # stack
    x_bins = np.arange(np.min(data[:,1]),np.max(data[:,1]+time_bin),time_bin)
    y_bins = np.arange(np.min(data[:,0]),np.max(data[:,0]+energy_bin),energy_bin)
    bin_data,x_edges,y_edges,_ = binned_statistic_2d(data[:,1], data[:,0], None, statistic='count', bins=[x_bins,y_bins])
    x_bins = np.arange(np.min(data[:,1]),np.max(data[:,1]),time_bin)



    for a in ranges: 
        low = a[0]
        high = a[1]
        
        xindmin = np.argmin(np.abs(y_bins-low))
        xindmax = np.argmin(np.abs(y_bins-high))

        print(xindmin, xindmax)
        #print(bin_data[(bin_data[:,0]>low)&(bin_data[:,0]<high)])
        #bin_data = np.sum(bin_data[(bin_data[:,1]>xindmin)&(bin_data[:,1]<xindmax)],axis=1)
        #bin_data = np.sum(bin_data[:,xindmin:xindmax],axis=1)
        bin_data2 = np.sum(bin_data[:,xindmin:xindmax],axis=1)
        norm = np.average(bin_data2[int(0.5*len(bin_data2)):int(0.75*len(bin_data2))])
        lab = (high+low)/2
        plt.plot(x_bins, bin_data2/norm,marker='.', label=str(lab), ls='none')
        plt.xlabel('time (ms)')
        plt.ylabel('Photon counts')
        plt.legend()



    plt.show()

plot_hslice([[1240,1250],[1200,1210],[835,845]])