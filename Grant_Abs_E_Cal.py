import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd

###Locating file path with proper dates and run numbers 
d = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
datdest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
today = "20221221"
rn = "0002"
#12/21      ##run0002       ###Declaring calibration states and science states 
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

for sta in scistates:
    histall = np.array(data.hist(np.arange(200, 9000, plotbinSize), "energy", states=sta))
    histall = histall.T
    stadat = pd.DataFrame(data=histall)
    stadat.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+'.csv', index=False)