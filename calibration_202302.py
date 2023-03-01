import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
#import ebit_util
import pandas as pd

#d = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
d = 'C:\\data\\TES_newcal'
#d = "C:\\data\\tesdata"
today = "20230216"
rn = "0000"
#datdest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022'
#datdest = 'C:\\data\\TES_Spectra_1eVbin'

savedat = True

#02/16 
calstates = ['A', 'B', 'C']

#02/17 run 0000
#calstates = ['A']

#run0001
# calstates = ['A', 'B', 'C']

fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)

data = ChannelGroup(fltoday)

defbinsize = 0.5
data.setDefaultBinsize(defbinsize)
data.learnResidualStdDevCut()
plotbinSize = 1


ds = data[1]
#ds.plotHist( np.arange(0, 60000, 10), "filtValue", coAddStates=False)   #states=None by default uses all states


ds.calibrationPlanInit("filtValue")

# #2/16 , run0000, chan1 
ds.calibrationPlanAddPoint(9790, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(11370, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(16660, "ClKAlpha", states=calstates)
ds.calibrationPlanAddPoint(20630, "KKAlpha", states=calstates)
ds.calibrationPlanAddPoint(27270, "TiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(37150, "FeKAlpha", states=calstates)


# #2/17 , run0000, chan1 
# ds.calibrationPlanAddPoint(9810, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11410, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16760, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20680, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27380, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37330, "FeKAlpha", states=calstates)


# #2/17 , run0001, chan1 
# ds.calibrationPlanAddPoint(9830, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11400, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16670, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20690, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27370, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37310, "FeKAlpha", states=calstates)


#ds.plotAvsB("relTimeSec", "filtValue", states=calstates)


data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
#data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
# Phase, drift, and time drift correct on pulses in the rough 400-2500 eV range
#data.cutAdd("Cut1", lambda energyRough: np.logical_and(energyRough > 400, energyRough < 3000), _rethrow=True)
data.learnPhaseCorrection(uncorrectedName="filtValue", correctedName = "filtValuePC", states=calstates)

data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValuePC", correctedName = "filtValuePCDC", states=calstates)
#data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC",)
#data.learnDriftCorrection(uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC",)

data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates, _rethrow=True) 
#data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates, _rethrow=True, cutRecipeName="cutForLearnDC",) 

data.calibrateFollowingPlan("filtValuePCDCTC", calibratedName="energy", overwriteRecipe=True)      



# external_trigger_filename =  os.path.join(d, f"{today}", f"{rn}", 
# f"{today}_run{rn}_external_trigger.bin")
# external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
# for ds in data.values():
#     ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)



#firstsci = scistates[0]
#firstsci = "F"
# energies = np.hstack([ds.getAttr("energy", firstsci) for ds in data.values()])
# seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]] for ds in data.values()])

#np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]][ds.getAttr("cutResidualStdDev", firstsci)] for ds in data.values()])

# for scistate in scistates: 
#     # energies = np.append(energies,np.hstack([ds.getAttr("energy", scistate) for ds in data.values()]))
#     # seconds_after_external_triggers = np.append(seconds_after_external_triggers,np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate)[0]] for ds in data.values()]))

#     energies = np.hstack([ds.getAttr("energy", scistate) for ds in data.values()])
#     seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate)[0]] for ds in data.values()])


#     dat2 = np.vstack((energies,seconds_after_external_triggers))
#     dat2 = dat2.T 
#     if savedat: 
#         np.save(datdest+'\\'+str(today)+'_'+str(rn)+'_'+str(scistate), dat2)

#data.plotHist(np.arange(0,10000,1),"energy", states=scistates, coAddStates=False)


# data.plotHist(np.arange(0,60000,10),"energy", states=scistates, coAddStates=True)
# ds.plotHist(np.arange(0,60000,10),"filtValue", states=scistates, coAddStates=True)
# data.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=calstates)
data.plotHist( np.arange(0,10000,1), "energy", coAddStates=True, states=calstates)
ds.diagnoseCalibration()
plt.show()


# for sta in scistates:
#     histall = np.array(data.hist(np.arange(500, 8000, plotbinSize), "energy", states=sta))
#     histall = histall.T
#     stadat = pd.DataFrame(data=histall)
    #stadat = stadat.transpose
    
    #stadat.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+'.csv', index=False)

    # energy = []
    # time = []
    # t_ext = []
    # for i in data:
    #     ds = data[i] 
    #     energy.extend(list(ds.getAttr('energy', sta)))
    #     time.extend(list(ds.getAttr('unixnano', sta)))
    #     t_ext.extend(list(ds.seconds_after_external_trigger[ds.getStatesIndicies(states=sta)[0]]))
    # #plist = np.array([energy, time, t_ext], dtype=object).T
    # plist = np.array([energy, time], dtype=object).T

    # #photlist = pd.DataFrame(data=plist, columns=['energy', 'time', 'external_trigger_time'])
    # photlist = pd.DataFrame(data=plist, columns=['energy', 'time'])
    
    #photlist.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+str('photonlist')+'.csv', index=False)