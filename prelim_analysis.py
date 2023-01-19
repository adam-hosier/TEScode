import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd

floc = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022')

date = str('202212')
day = str('15')
runnum = str('0001')
statelist = ['G', 'K', 'O', 'Q', 'U', 'Y', 'Z', 'AB']
coAdd = False
dfall = dict()
minenergy = 500
maxenergy = 5000
binsize = 1
numbins = int(np.round((maxenergy-minenergy)/binsize))


for s in statelist: 
    state = str(s)
    dfall[state] = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')
    df = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')

    counts, bin_edges = np.histogram(dfall[state]['energy'], bins=numbins, range=(minenergy, maxenergy))
    dfall[state+str(' counts')]= counts
    dfall[state+str(' bin_edges')] = bin_edges

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

#plt.figure()
plt.ylabel('Counts per '+str(binsize)+' eV bin')
#plt.ylim(bottom=0, top=1.1*np.max(counts))
#plt.xlim(left=3075, right=3175)
plt.xlabel('energy (eV)')
plt.title(date+day+'_'+runnum)
#plt.title(date+day+'_'+runnum+'_'+state)
for state in statelist: 
    plt.plot(dfall[state+str(' bin_edges')][:-1], dfall[state+str(' counts')], label=state)
plt.legend()
plt.show() 


