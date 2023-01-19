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
state = str('G')

df = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')

energy = df['energy']
time = df['time']

minenergy = 500
maxenergy = 5000
binsize = 1
numbins = int(np.round((maxenergy-minenergy)/binsize))

counts, bin_edges = np.histogram(energy, bins=numbins, range=(minenergy, maxenergy))


##actual histogram
# plt.figure()
# plt.ylabel('Counts per '+str(binsize)+' eV bin')
# plt.ylim(bottom=0, top=1.1*np.max(counts))
# #plt.xlim(left=3075, right=3175)
# plt.xlabel('energy (eV)')
# plt.title(date+day+'_'+runnum+'_'+state)
# plt.hist(energy, bins=numbins, range=(minenergy, maxenergy))
# plt.show() 

plt.figure()
plt.ylabel('Counts per '+str(binsize)+' eV bin')
plt.ylim(bottom=0, top=1.1*np.max(counts))
plt.xlim(left=3075, right=3175)
plt.xlabel('energy (eV)')
plt.title(date+day+'_'+runnum+'_'+state)
plt.plot(bin_edges[:-1], counts, label=state)
plt.legend()
plt.show() 


