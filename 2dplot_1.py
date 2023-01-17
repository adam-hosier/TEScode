import math
import time
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import ebit_util

cook_time = 600.95 # cook time in ms
stack_period = 62 # stack period (duty cycle) in ms

v_ramp_t = [0,3,13,15,25,28,62] # DT timing points
v_ramp_v = [5000,3200,2600,2600,3200,5000,5000] # corresponding DT voltages in V

plot_energy_bounds = [0,10000] # [upper, lower] photon energy axis limits (eV)
plot_voltage_bounds = [2200,4990] # [upper, lower] beam energy axis limits (eV)

energy_bin = 2 # photon energy bins size (eV)
voltage_bin = 5 # beam energy bins size (eV)

I = 50 # ebeam current in mA

count_to_print = [[2220,3100,'He'],[2220,3200,'H']] # array of x,y values to count and label
count_radius = 20 # radius of count_to_print to be integrated

x_lims = [1000,5000]
y_lims = [1000,5000]

d = "/home/pcuser/data"
date = "20221215"                       #date
rn = "0001"                             #run number, can be found in dastard command 
state1 = "B"
fl = getOffFileListFromOneFile(os.path.join(d, f"{date}", f"{rn}", 
f"{date}_run{rn}_chan1.off"), maxChans=300)
data = ChannelGroup(fl)             #all of the channels 
data.setDefaultBinsize(0.4)

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

#calibrate and inspect one channel
ds = data[1]                       #selection of individual channel
ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=state1 )   #states=None by default uses all states


ds.learnDriftCorrection(overwriteRecipe=True, states=state1)
#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
ds.calibrationPlanInit("filtValueDC")
ds.calibrationPlanAddPoint(19431, "KKAlpha", states=state1)
ds.calibrationPlanAddPoint(35592, "FeKAlpha", states=state1)
# ds.calibrationPlanAddPoint(37150, "FeKBeta", states=state1)
ds.calibrationPlanAddPoint(25950, "TiKAlpha", states=state1)
ds.calibrationPlanAddPoint(8947, "AlKAlpha", states=state1)
ds.calibrationPlanAddPoint(15549, "ClKAlpha", states=state1)
ds.calibrationPlanAddPoint(10447, "SiKAlpha", states=state1)


ds.plotHist( np.arange(0, 10000, 1), "energyRough", coAddStates=False, states=state1 )   #states=None by default uses all states
#ds.plotAvsB("relTimeSec", "pretriggerMean")
#ds.plotAvsB("relTimeSec", "energyRough")
#plt.ylim(1100,1600)
ds.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)


#ds.plotHist( np.arange(0,8500,2), "energyRough", coAddStates=False, states=state1)

ds.calibrateFollowingPlan("filtValuePC")
ds.diagnoseCalibration()
ds.plotAvsB("filtPhase", "energy")
ds.plotAvsB("pretriggerMean","energy")

data.learnDriftCorrection(overwriteRecipe=True, states=state1)
data.alignToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
data.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
data.calibrateFollowingPlan("filtValuePC")

def refresh_and_print_if_worked():
    old = np.sum([len(dslocal) for dslocal in data.values()])
    data.refreshFromFiles()
    new = np.sum([len(dslocal) for dslocal in data.values()])
    if new > old:
        print(f"refresh worked, now have {new} pulses vs {old} before")
    else:
        print(f"refresh FAILED, now have {new} pulses vs {old} before")

#realtime 2d specta, re-run this cell to refresh
external_trigger_filename =  os.path.join(d, f"{date}",f"{rn}", f"{date}_run{rn}_external_trigger.bin")
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)

states_for_2dhist = ["G", "K", "O", "Q"]
energies = np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()])

print(data.shortName)


def plot2d(energies, seconds_after_external_triggers):

    data = np.vstack((energies,seconds_after_external_triggers))
    data = data.T
    np.save('/home/pcuser/Desktop/TES-GUI/20221215_0001_GKOQ', data)
    data = data[(data[:,0] > plot_energy_bounds[0]) & (data[:,0] <= plot_energy_bounds[1])] # truncate data out of energy bounds
    data[:,1] = (data[:,1] - cook_time) % stack_period # stack
    voltage_func = interp1d(v_ramp_t,v_ramp_v,kind='linear') # interpolated DT voltage function
    data = np.column_stack((data, data[:,1])) #adding a new column that will become DT voltage
    data[:,2] = voltage_func(data[:,2]) # voltage to time map
    data = data[(data[:,2] > plot_voltage_bounds[0]) & (data[:,2] <= plot_voltage_bounds[1])] # truncate data out of voltage bounds
    data = sc_correct(data,I) # apply space charge correction
    x_bins = np.arange(np.min(data[:,2]),np.max(data[:,2]+voltage_bin),voltage_bin)
    y_bins = np.arange(np.min(data[:,0]),np.max(data[:,0]+energy_bin),energy_bin)
    bin_data,x_edges,y_edges,_ = binned_statistic_2d(data[:,2], data[:,0], None, statistic='count', bins=[x_bins,y_bins])
    bin_data = np.ma.masked_array(bin_data,bin_data<1)
    fig = plt.figure()
    plt.xlabel('Beam Energy (eV)')
    plt.ylabel('Photon Energy (eV)')
    for count_loc in count_to_print:
        print(f'{count_loc[2]} {counts(count_loc[0], count_loc[1], count_radius, x_edges, y_edges, bin_data)} counts')
    plt.pcolormesh(x_edges,y_edges,bin_data.T)
    #plt.xlim(x_lims[0],x_lims[1])
    #plt.ylim(y_lims[0],y_lims[1])
    fig.canvas.draw()
    plt.colorbar()
    
plot2d(energies, seconds_after_external_triggers)