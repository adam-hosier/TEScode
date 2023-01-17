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

v_ramp_t = [0,3,13,15,25,28,62] # DT timing points
v_ramp_v = [5000,3200,2600,2600,3200,5000,5000] # corresponding DT voltages in V

plot_energy_bounds = [0,10000] # [upper, lower] photon energy axis limits (eV)
plot_voltage_bounds = [2600,4940] # [upper, lower] beam energy axis limits (eV)

energy_bin = 2 # photon energy bins size (eV)
time_bin = 5 # time bins size (ms)


def plot_hslice(low,high):

    data = np.load('/home/pcuser/Desktop/TES-GUI/20221215_0001_9states.npy')
    
    data = data[(data[:,0] > plot_energy_bounds[0]) & (data[:,0] <= plot_energy_bounds[1])] # truncate data out of energy bounds
    data[:,1] *= 1000
    data = data[data[:,1] < stack_period*cycles]
    data[:,1] = data[:,1] % stack_period # stack

    x_bins = np.arange(np.min(data[:,1]),np.max(data[:,1]+time_bin),time_bin)
    y_bins = np.arange(np.min(data[:,0]),np.max(data[:,0]+energy_bin),energy_bin)
    bin_data,x_edges,y_edges,_ = binned_statistic_2d(data[:,1], data[:,0], None, statistic='count', bins=[x_bins,y_bins])

    bin_data = np.sum(bin_data[(bin_data[:,0]>low)&(bin_data[:,0]<high)],axis=1)
    plt.plot(x_bins,bin_data,marker='.')
    plt.show()

plot_hslice(3000,3200)