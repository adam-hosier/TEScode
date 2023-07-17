import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd


d = 'C:\\data\\TimingData'
SoI = ['K', 'M'] #states of interest for slicing 
files = dict() 

for states in SoI: 
    path1 = os.path.join(d, "20221221_0002_"+str(states)+".npy")
    files[str(states)] = np.load(path1)

plotbinSize = 1
#np.arange(200, 9000, plotbinSize)

minE = 200 
maxE = 9000
numbin = int((maxE-minE)/plotbinSize)
counts, bins = np.histogram(files['K'], bins=numbin, range=(minE, maxE))


# plt.stairs(counts, bins, fill=False)
# plt.show()
timebinsize = 1/1000         #bin size for timing in s 
timebin = np.arange(0, 3, timebinsize)
## Co-like 1182 eV \\ Ni-like 1144
intwidth = 5

nilikecounts = []
niliketime = []
colikecounts = []
coliketime = []

nicounts = []
nitime = []
cocounts = []
cotime = []

files['K'] = np.delete(files['K'], np.where(files['K']<0)[0], axis=0)


hist, xedges, yedges = np.histogram2d(files['K'][:,1], files['K'][:,0], bins=((int((3-0)/timebinsize)), 500))
print('done 2dhist')

plt.hist2d(files['K'][:,1], files['K'][:,0], bins=int((3-0)/timebinsize))
plt.colorbar()
plt.show()
# for states in SoI: 

#     cocounter = 0
#     nicounter = 0
#     for i in range(np.shape(files[str(states)])[0]):

#         if files[str(states)][i,0] <= 1182+intwidth and files[str(states)][i,0] >= 1182-intwidth:
#             colikecounts.append(files[str(states)][i,0])
#             coliketime.append(files[str(states)][i,1])
#         elif files[str(states)][i,0] <= 1144+intwidth and files[str(states)][i,0] >= 1144-intwidth:
#             nilikecounts.append(files[str(states)][i,0])
#             niliketime.append(files[str(states)][i,1])
