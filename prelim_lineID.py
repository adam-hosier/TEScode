import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
#import ebit_util
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import binned_statistic_2d
from scipy import stats
from fit_utils import MultiPeakGaussian

floc = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022')
#floc = str('C:\\data\\calibratedTES_Dec2022')
ddest = str('C:\\data\\Line_ID_Nd')
date = str('202212')
day = str('21')
runnum = str('0002')
statelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']

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

energy = df['energy']
time = df['time']

# counts, bin_edges = np.histogram(energy, bins=numbins, range=(minenergy, maxenergy))

###################################
# plt.figure()
# plt.ylabel('Counts per '+str(binsize)+' eV bin')
# #plt.ylim(bottom=0, top=1.1*np.max(counts))
# #plt.xlim(left=3075, right=3175)
# plt.xlabel('energy (eV)')
# plt.title(date+day+'_'+runnum)
# #plt.title(date+day+'_'+runnum+'_'+state)
# for state in statelist: 
#     plt.plot(dfall[state+str(' bin_edges')][:-1], dfall[state+str(' counts')], label=state)
# plt.legend()
# #plt.show() 
# ##################################
state = 'T'
arry = dfall[state+str(' counts')]
arrx = dfall[state+str(' bin_edges')][:-1]
res_dict = dict()
a = MultiPeakGaussian(arr = arry, xs = arrx, num_peaks=60, num_poly=3)
a.fit(return_dict = res_dict, same_sigma=True, function='voigt')
t1 = res_dict['rez']
a.plot_fit(normalize_background=True)

b = t1[1]
c = np.zeros([1,9])
for i in range(np.shape(t1)[0]):

    tempa = np.zeros([1,1])
    for l in range(np.shape(t1)[1]):
        temp = np.array([t1[i][l]])

        tempa = np.hstack((tempa,temp))
    
    c = np.vstack((c, tempa))

c = c[1:,:]
c = c[:,1:]

lineIDdata = pd.DataFrame(data=c, columns=['center [eV]', 'center unc [eV]','height [counts]','height unc [counts]', 
    'std dev [eV]','std dev unc [eV]', 'FWHM [eV]', 'FWHM unc [eV]'])

lineIDdata.to_csv(ddest+'/'+str(date)+str(day)+'_'+str(runnum)+'_'+str(state)+'.csv', index=True)