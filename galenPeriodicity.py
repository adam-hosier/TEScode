#import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd
from scipy.fft import fft, fftfreq

d = os.path.dirname(os.path.abspath(__file__))
today = "20221220"
rn = "0000"
datdest = os.path.dirname(os.path.abspath(__file__))
#12/20 RUN_0000
calstates = ["A", "B", "K"]    
fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan4.off"), maxChans=300)
external_trigger_filename =  os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_external_trigger.bin")

data = ChannelGroup(fltoday)
data.setDefaultBinsize(0.5)
ds = data[0]
ds.calibrationPlanInit("filtValue")
#12/20
# ##run0000       chan 4 filt values          
ds.calibrationPlanAddPoint(9819, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(11400, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(16702, "ClKAlpha", states=calstates)
ds.calibrationPlanAddPoint(20733, "KKAlpha", states=calstates)
data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
data.learnPhaseCorrection(uncorrectedName="filtValue", correctedName = "filtValuePC", states=calstates)
data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValuePC", correctedName = "filtValuePCDC", states=calstates)
data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates) 
data.calibrateFollowingPlan("filtValuePCDCTC", calibratedName="energy", overwriteRecipe=True)      


external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)

firstsci = calstates[0]
energies = np.hstack([ds.getAttr("energy", firstsci) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]] for ds in data.values()])


for scistate in calstates: 
    energies = np.hstack([ds.getAttr("energy", scistate) for ds in data.values()])
    seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate)[0]] for ds in data.values()])
    dat2 = np.vstack((energies,seconds_after_external_triggers))
    dat2 = dat2.T 
    #np.save(datdest+'\\'+str(today)+'_'+str(rn)+'_'+str(scistate), dat2)
    np.save(os.path.join(d, f"{today}", f"{rn}", f"{scistate}"), dat2)

plt.close() 
##Calibration complete, moving onto plotting FFT
plt.rcParams.update({'font.size': 16})
dir = datdest
file = '20221220_0000_'
states = ['A','B','K']
binsize = .001
e_bounds = [2610,2630]
t_bounds = [0,1]
data = np.empty((1,2))
for state in states:
    data = np.vstack((data,np.load(os.path.join(dir, f"{today}", f"{rn}", f"{state}.npy"))))
data = data[(data[:,0]>e_bounds[0]) & (data[:,0]<e_bounds[1])]
data = data[(data[:,1]>t_bounds[0]) & (data[:,1]<t_bounds[1])]
bin_edges = np.arange(t_bounds[0], t_bounds[1]+binsize, binsize)
counts, _ = np.histogram(data[:,1], bins=bin_edges)
bin_centers = bin_edges[:-1]+binsize/2
N = len(counts)
T = binsize
xf = fftfreq(N,T)[1:N//2]
plt.plot(xf,2/N*np.abs(fft(counts)[1:N//2]))
plt.xlabel('Freq [Hz]')
plt.ylabel('Magnitude (arb)')
plt.minorticks_on()
plt.show()