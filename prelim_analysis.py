import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import binned_statistic_2d
from scipy import stats
# floc = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022')

# date = str('202212')
# day = str('15')
# runnum = str('0001')
# statelist = ['G', 'K', 'O', 'Q', 'U', 'Y', 'Z', 'AB']
# coAdd = False
# dfall = dict()
# minenergy = 500
# maxenergy = 5000
# binsize = 1
# numbins = int(np.round((maxenergy-minenergy)/binsize))


# for s in statelist: 
#     state = str(s)
#     dfall[state] = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')
#     df = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')

#     counts, bin_edges = np.histogram(dfall[state]['energy'], bins=numbins, range=(minenergy, maxenergy))
#     dfall[state+str(' counts')]= counts
#     dfall[state+str(' bin_edges')] = bin_edges

# energy = df['energy']
# time = df['time']


# counts, bin_edges = np.histogram(energy, bins=numbins, range=(minenergy, maxenergy))


##actual histogram
# plt.figure()
# plt.ylabel('Counts per '+str(binsize)+' eV bin')
# plt.ylim(bottom=0, top=1.1*np.max(counts))
# #plt.xlim(left=3075, right=3175)
# plt.xlabel('energy (eV)')
# plt.title(date+day+'_'+runnum+'_'+state)
# plt.hist(energy, bins=numbins, range=(minenergy, maxenergy))
# plt.show() 



###################################
# #plt.figure()
# plt.ylabel('Counts per '+str(binsize)+' eV bin')
# #plt.ylim(bottom=0, top=1.1*np.max(counts))
# #plt.xlim(left=3075, right=3175)
# plt.xlabel('energy (eV)')
# plt.title(date+day+'_'+runnum)
# #plt.title(date+day+'_'+runnum+'_'+state)
# for state in statelist: 
#     plt.plot(dfall[state+str(' bin_edges')][:-1], dfall[state+str(' counts')], label=state)
# plt.legend()
# plt.show() 
# ##################################

dat2d_F_loc = 'C:\\data\\calibrated_data\\20221220_0001_F.npy'
dat2d_N_loc = 'C:\\data\\calibrated_data\\20221220_0001_N.npy'
dat2d_F = np.load(dat2d_F_loc)
dat2d_N = np.load(dat2d_N_loc)
energies_F = dat2d_F[:,0]
seconds_after_external_triggers_F = dat2d_F[:,1]
energies_N = dat2d_N[:,0]
seconds_after_external_triggers_N = dat2d_N[:,1]


energies = np.append(energies_F, energies_N)
seconds_after_external_triggers = np.append(seconds_after_external_triggers_F, seconds_after_external_triggers_N)

energies = energies_N
seconds_after_external_triggers = seconds_after_external_triggers_N
plt.figure()
plt.hist2d(seconds_after_external_triggers, 
    energies, bins=(np.arange(0,3.5,0.01), np.arange(500,2500,1)),
    cmin = 1, cmax = 100)
plt.xlabel("time since external trigger (s)")
plt.ylabel("energy(eV)")
#plt.colormesh(norm=mcolors.PowerNorm(0.5))
#plt.title(f"{data.shortName}, states={states_for_2dhist}")
plt.colorbar()


time_bin = 0.005
energy_bin = 5
low = 800
high = 850
#def plot_hslice(low,high):

data2 = np.load(dat2d_N_loc)

x_bins = np.arange(np.min(data2[:,1]),np.max(data2[:,1]+time_bin),time_bin)
#y_bins = np.arange(np.min(data2[:,0]),np.max(data2[:,0]+energy_bin),energy_bin)
y_bins = np.arange(0,np.max(data2[:,0]+energy_bin),energy_bin)

xindmin = np.argmin(np.abs(y_bins-low))
xindmax = np.argmin(np.abs(y_bins+high))

bin_data,x_edges,y_edges,_ = binned_statistic_2d(data2[:,1], data2[:,0], None, statistic='count', bins=[x_bins,y_bins])

bin_data = np.sum(bin_data[(bin_data[:,0]>low)&(bin_data[:,0]<high)],axis=0)
bin_data2 = np.sum(bin_data[xindmin:xindmax], axis=0)


ne = []
nt = []
for x in range(len(energies)):
    if energies[x] <= high and energies[x] > low:
        ne.append(energies[x])
        nt.append(seconds_after_external_triggers[x])
    else: 
        pass 

y_bins2 = np.arange(low, high+energy_bin, energy_bin)
# ne = np.array(ne)
# nt = np.array(nt)
bd2, xedg, yedg = stats.binned_statistic(nt, ne, statistic='count', bins=100)



# print(bin_data2)
# print(x_edges)
# print(y_edges)


plt.plot(xedg[:-1]+time_bin/2,bd2,marker='.', ls='none')

#plt.plot(x_bins,bin_data,marker='.', ls='none')
# plt.plot(x_bins, bin_data2, marker='.', ls='none')
# plt.show()

#plot_hslice(800,850)

