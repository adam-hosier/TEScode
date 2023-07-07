import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd
#plt.ion()
d = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
today = "20221221"
rn = "0002"
datdest = 'C:\\data\\TimingData'

#12/21      ##run0002
calstates = ["A", "I", "O", "AH"]
scistates = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]
fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)
data = ChannelGroup(fltoday)
defbinsize = 0.5
data.setDefaultBinsize(defbinsize)
plotbinSize = 1
ds = data[1]

ds.calibrationPlanInit("filtValue")
#12/21              
# ##run0002             #chan 1
ds.calibrationPlanAddPoint(8190, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(9550, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(14070, "ClKAlpha", states=calstates)
#ds.calibrationPlanAddPoint(15080, "ClKBeta", states=calstates)
ds.calibrationPlanAddPoint(17570, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(23279, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(31720, "FeKAlpha", states=calstates)

data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
data.learnPhaseCorrection(uncorrectedName="filtValue", correctedName = "filtValuePC", states=calstates)
data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValuePC", correctedName = "filtValuePCDC", states=calstates)
data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates) 
data.calibrateFollowingPlan("filtValuePCDCTC", calibratedName="energy", overwriteRecipe=True)      

external_trigger_filename =  os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_external_trigger.bin")
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)

firstsci = scistates[0]
energies = np.hstack([ds.getAttr("energy", firstsci) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]] for ds in data.values()])

for scistate in scistates: 
    energies = np.hstack([ds.getAttr("energy", scistate) for ds in data.values()])
    seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate)[0]] for ds in data.values()])
    dat2 = np.vstack((energies,seconds_after_external_triggers))
    dat2 = dat2.T 
    np.save(datdest+'\\'+str(today)+'_'+str(rn)+'_'+str(scistate), dat2)