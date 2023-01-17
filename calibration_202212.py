import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd

d = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
today = "20221220"
rn = "0001"
datdest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022'


#12/14 
#calstates = ["B", "C"]

#12/15 
# calstates = ["A", "B", "I", "M", "S", "W", "AF"]
# scistates = ["G", "K", "O", "Q", "U" ,"Y", "Z", "AB", "AD", "AH", "AO"]

#12/16
# calstates = ["A", "B", "F", "J", "N", "R"]    ##run0000
# scistates = ["E", "H", "L", "P"]

# calstates = ["B", "F"]   ##run0001
# scistates = ["D"]

#12/19          ## run 0000
# calstates = ["A", "B", "C", "D", "AC"]
# scistates = ["H", "K", "P", "W", "Y", "AA", "R", "T", "U"]


#12/20 
# calstates = ["A", "B", "K"]       ##run0000
# scistates = ["H", "J", "L", "P"]

calstates = ["B", "D", "J"]        ##run0001 
scistates = ["F", "N"]
scistates = ["F"]

#12/21
#calstates = ["A", "B", "I", "O", "AH"]
#scistates = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]



fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)


data = ChannelGroup(fltoday)

defbinsize = 0.2
data.setDefaultBinsize(defbinsize)



ds = data[1]
#ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=calstates)   #states=None by default uses all states


ds.calibrationPlanInit("filtValue")

#12/15 
# ds.calibrationPlanAddPoint(8950, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(15560, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(19440, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(25950, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(35570, "FeKAlpha", states=calstates)


# #12/16   
# #run0000
# ds.calibrationPlanAddPoint(8660, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10049, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(14785, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(18292, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(24085, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(32678, "FeKAlpha", states=calstates)


# #12/16  
# #run0001
# ds.calibrationPlanAddPoint(9820, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11440, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16828, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20848, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27494, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37210, "FeKAlpha", states=calstates)



# #12/19 
# ds.calibrationPlanAddPoint(8970, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(15377, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(19030, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(25190, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(34610, "FeKAlpha", states=calstates)


#12/20
##run0000
# ds.calibrationPlanAddPoint(8994, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10454, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(15398, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(19020, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(25230, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(34600, "FeKAlpha", states=calstates)

##run0001
ds.calibrationPlanAddPoint(8970, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(10450, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(15390, "ClKAlpha", states=calstates)
ds.calibrationPlanAddPoint(19020, "KKAlpha", states=calstates)
ds.calibrationPlanAddPoint(25220, "TiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(34610, "FeKAlpha", states=calstates)



#12/21              
# ##run0002
# ds.calibrationPlanAddPoint(8190, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(9550, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(14090, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(17550, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(23270, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(31720, "FeKAlpha", states=calstates)
#ds.plotAvsB("relTimeSec", "filtValue", states=calstates)

data.learnResidualStdDevCut()
data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
#data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
# Phase, drift, and time drift correct on pulses in the rough 400-2500 eV range
#data.cutAdd("Cut1", lambda energyRough: np.logical_and(energyRough > 400, energyRough < 3000), _rethrow=True)
data.learnPhaseCorrection(uncorrectedName="filtValue", correctedName = "filtValuePC", states=calstates)

#data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValuePC", correctedName = "filtValuePCDC", states=calstates, cutRecipeName="cutForLearnDC")
#data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC",)
#data.learnDriftCorrection(uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC",)

#data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC", _rethrow=True) 
#data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates, _rethrow=True, cutRecipeName="cutForLearnDC",) 

data.calibrateFollowingPlan("filtValuePC", calibratedName="energy", overwriteRecipe=True)               #         
#data.plotHist(np.arange(0,10000,1),"energy", states=scistates, coAddStates=False)


# data.plotHist(np.arange(0,60000,10),"energy", states=scistates, coAddStates=True)
# ds.plotHist(np.arange(0,60000,10),"filtValue", states=scistates, coAddStates=True)

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