import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd

d = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
today = "20221221"
rn = "0002"

#12/14 
#calstates = ["B", "C"]

#12/15 
# calstates = ["A", "B", "I", "M", "S", "W", "AF"]

#12/16
#calstates = ["A", "B", "F", "J", "N", "R"]

#12/19
#calstates = ["A", "B", "C", "D", "AC"]

#12/20 0000 
#calstates = ["A", "B", "K", "N"]

#0001 
#calstates = ["B", "D", "J"]

#12/21 
#calstates = ["A", "B", "I", "O", "AH"]


#12/21
calstates = ["A", "B", "I", "O", "AH"]


datdest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022'

fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)


data = ChannelGroup(fltoday)

defbinsize = 0.5
data.setDefaultBinsize(defbinsize)



ds = data[1]
#ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=calstates)   #states=None by default uses all states
# plt.show()
# asdfas

#ds.learnTimeDriftCorrection(states=calstates
#)
ds.calibrationPlanInit("filtValue")
ds.learnDriftCorrection(overwriteRecipe=True, states=calstates)
#ds.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states


#12/15 
# ds.calibrationPlanAddPoint(8950, "AlKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(15560, "ClKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(19440, "KKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(25950, "TiKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(35570, "FeKAlpha", states=calstates
#)


# #12/16 
# ds.calibrationPlanAddPoint(8660, "AlKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(10049, "SiKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(14785, "ClKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(18292, "KKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(24085, "TiKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(32678, "FeKAlpha", states=calstates
#)


# #12/19 
# ds.calibrationPlanAddPoint(8970, "AlKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(15377, "ClKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(19030, "KKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(25190, "TiKAlpha", states=calstates
#)
# ds.calibrationPlanAddPoint(34610, "FeKAlpha", states=calstates
#)


#12/21
ds.calibrationPlanAddPoint(8190, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(9550, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(14090, "ClKAlpha", states=calstates)
ds.calibrationPlanAddPoint(17550, "KKAlpha", states=calstates)
ds.calibrationPlanAddPoint(23270, "TiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(31720, "FeKAlpha", states=calstates)

scistates = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]

data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,35000,10), attr="filtValue", states=calstates)
# Phase, drift, and time drift correct on pulses in the rough 400-2500 eV range
data.cutAdd("cutForLearnDC", lambda energyRough: np.logical_and(energyRough > 400, energyRough < 2500), setDefault=False, _rethrow=True)
data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC", states=calstates)
data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValuePC", correctedName = "filtValuePCDC", states=calstates, cutRecipeName="cutForLearnDC")
data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates,
                            cutRecipeName="cutForLearnDC", _rethrow=True) 

data.calibrateFollowingPlan("filtValuePCDCTC", calibratedName="energy", overwriteRecipe=True)                        
data.plotHist(np.arange(0,10000,1),"energy", states=scistates)



#ds.learnPhaseCorrection(uncorrectedName = "filtValueDC", linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
#ds.calibrateFollowingPlan("filtValuePC")
#ds.learnPhaseCorrection(states=calstates
#)

#ds.calibrateFollowingPlan("filtValueDC")

#data.learnDriftCorrection(overwriteRecipe=True, states=calstates)

#data.aligniiToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
#data.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
#data.learnPhaseCorrection(states=calstates
#)
#data.calibrateFollowingPlan("filtValueDC")



#data.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=scistates)

# for sta in scistates:
#     histall = np.array(data.hist(np.arange(0, 14000, defbinsize), "energy", states=sta))
#     stadat = pd.DataFrame(data=histall)
#     #stadat = stadat.transpose
#     stadat.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+'.csv', index=False)

#     energy = []
#     time = []
#     for i in data:
#         ds = data[i] 
#         energy.extend(list(ds.getAttr('energy', sta)))
#         time.extend(list(ds.getAttr('unixnano', sta)))
#     plist = np.array([energy, time], dtype=object).T

#     photlist = pd.DataFrame(data=plist, columns=['energy', 'time'])
#     photlist.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+str('photonlist')+'.csv', index=False)

