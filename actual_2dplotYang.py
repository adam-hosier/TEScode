import math
import time
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

cook_time = 1000.95 # cook time in ms
stack_period = 60 # stack period (duty cycle) in ms
cycles = 33

v_ramp_t = [0,3,13,15,25,28,60] # DT timing points
v_ramp_v = [9500,6900,6600,6600,6900,9500,9500] # corresponding DT voltages in V

plot_energy_bounds = [6000,8000] # [upper, lower] photon energy axis limits (eV)
plot_voltage_bounds = [6520,9000] # [upper, lower] beam energy axis limits (eV)

energy_bin = 20 # photon energy bins size (eV)
voltage_bin = 5 # beam energy bins size (eV)

I = 105 # ebeam current in mA

count_to_print = [[6697,6694,'He'],[2220,3200,'H']] # array of x,y values to count and label
count_radius = 20 # radius of count_to_print to be integrated


fload = '/home/pcuser/Desktop/TES-GUI/Yang/20221216_0000_H.npy'
fload2 = '/home/pcuser/Desktop/TES-GUI/Yang/20221216_0001_D.npy'
def sc_correct(vdata, e_curr): # data_arr and ebeam current in mA
    I = e_curr * 10**(-3) # beam current
    k = 8.99 * 10**(9) # Coulomb constant
    m = 9.11 * 10**(-31) # e- mass
    r = 2.5 * 10**(-3) # DT radius
    re = 200 * 10**(-6) # beam radius
    j2ev = 1.602 * 10**(-19) # J/eV aka e- charge
    f = 1 # multiplicative fudge factor, this should include the neutralization factor

    const = I*k/(math.sqrt(2/m))*(1+2*math.log(r/re)) *f

    vdata[:,2] -= const * (j2ev * vdata[:,2])**(-1/2)

    return vdata
def counts(x,y,r,xbins,ybins,vals):
    xshift = (xbins[1]-xbins[0])/2 # shift the bin values to the center
    yshift = (ybins[1]-ybins[0])/2
    xbins_shift = xbins + xshift
    ybins_shift = ybins + yshift

    xmatch = xbins_shift[(xbins_shift < x+r) & (xbins_shift > x-r)]
    ymatch = ybins_shift[(ybins_shift < y+r) & (ybins_shift > y-r)]
    pairs = [item for sublist in [[[i,j] for i in xmatch] for j in ymatch] for item in sublist]
    pairmatch = [i for i in pairs if ((i[0]-x)**2 + (i[1]-y)**2)**(.5) < r]

    pairmatch = [[i[0]-xshift,i[1]-yshift] for i in pairmatch] # shift back to bottom left
    indexpairs = [[np.argmin(abs(i[0]-xbins)),np.argmin(abs(i[1]-ybins))] for i in pairmatch]

    return np.nansum([vals[i[0],i[1]] for i in indexpairs])

def plot2d():

    data = np.load(fload)
    data2 = np.load(fload2)
    data = np.vstack((data,data2))
    data = data[(data[:,0] > plot_energy_bounds[0]) & (data[:,0] <= plot_energy_bounds[1])] # truncate data out of energy bounds
    data[:,1] *= 1000
    data = data[data[:,1] < stack_period*cycles]
    data[:,1] = data[:,1] % stack_period # stack
    voltage_func = interp1d(v_ramp_t,v_ramp_v,kind='linear') # interpolated DT voltage function
    x = np.linspace(0,62,100)
    #plt.scatter(data[:,1], data[:,0],marker='.',s=1)
    data = np.column_stack((data, data[:,1])) #adding a new column that will become DT voltage
    data[:,2] = voltage_func(data[:,2]) # voltage to time map

    #data = sc_correct(data,I) # apply space charge correction
    data = data[(data[:,2] > plot_voltage_bounds[0]) & (data[:,2] <= plot_voltage_bounds[1])] # truncate data out of voltage bounds
    
    x_bins = np.arange(np.min(data[:,2]),np.max(data[:,2]+voltage_bin),voltage_bin)
    y_bins = np.arange(np.min(data[:,0]),np.max(data[:,0]+energy_bin),energy_bin)
    bin_data,x_edges,y_edges,_ = binned_statistic_2d(data[:,2], data[:,0], None, statistic='count', bins=[x_bins,y_bins])
    bin_data = np.ma.masked_array(bin_data,bin_data<1)
    fig = plt.figure()
    plt.title('20221216_run0001 105 chans, states=D')
    plt.xlabel('Beam Energy (eV)')
    plt.ylabel('Photon Energy (eV)')
    for count_loc in count_to_print:
        print(f'{count_loc[2]} {counts(count_loc[0], count_loc[1], count_radius, x_edges, y_edges, bin_data)} counts')
    plt.pcolormesh(x_edges,y_edges,bin_data.T, norm=mcolors.PowerNorm(0.5))
    plt.colorbar()
    plt.show()
    

def slices(x_range, step):
    data = np.load(fload)
    data2 = np.load(fload2)
    data = np.vstack((data,data2))

    data = data[(data[:,0] > plot_energy_bounds[0]) & (data[:,0] <= plot_energy_bounds[1])] # truncate data out of energy bounds
    data[:,1] *= 1000
    data = data[data[:,1] < stack_period*cycles]
    data[:,1] = data[:,1] % stack_period # stack
    voltage_func = interp1d(v_ramp_t,v_ramp_v,kind='linear') # interpolated DT voltage function
    x = np.linspace(0,62,100)
    #plt.scatter(data[:,1], data[:,0],marker='.',s=1)
    data = np.column_stack((data, data[:,1])) #adding a new column that will become DT voltage
    data[:,2] = voltage_func(data[:,2]) # voltage to time map

    #data = sc_correct(data,I) # apply space charge correction
    data = data[(data[:,2] > plot_voltage_bounds[0]) & (data[:,2] <= plot_voltage_bounds[1])] # truncate data out of voltage bounds
    
    x_bins = np.arange(np.min(data[:,2]),np.max(data[:,2]+voltage_bin),voltage_bin)
    y_bins = np.arange(np.min(data[:,0]),np.max(data[:,0]+energy_bin),energy_bin)
    bin_data,x_edges,y_edges,_ = binned_statistic_2d(data[:,2], data[:,0], None, statistic='count', bins=[x_bins,y_bins])

    fig,ax = plt.subplots()
    y= np.arange(np.min(data[:,0]),np.max(data[:,0]),energy_bin)
    slicess = np.arange(x_range[0],x_range[1],step)
    step = slicess[1]-slicess[0]
    colors = np.linspace(0,1,len(slicess))
    for i,color in zip(slicess,colors):
        slice_range_low = np.argmin(abs(x_bins-(i)))
        slice_range_high =  np.argmin(abs(x_bins-(i+step)))
        ax.plot(y,np.sum(bin_data[slice_range_low:slice_range_high,:],axis=0),color=plt.cm.plasma(color),label=i)
    plt.legend()
    plt.show()

plot2d()
#slices([7000,9000],500)