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

#12/15 ##run0001        
# calstates = ["A", "B", "I", "M", "S", "W", "AF"]
# scistates = ["G", "K", "O", "Q", "U" ,"Y", "Z", "AB", "AD", "AH", "AO"]

#12/16
# calstates = ["A", "B", "F", "J", "N", "R"]    ##run0000
# scistates = ["E", "H", "L", "P"]

# calstates = ["B", "F"]   ##run0001
# #calstates = ["F"]

# scistates = ["D"]

#12/19          ## run 0000
# calstates = ["A", "B", "C", "D", "AC"]
# scistates = ["H", "K", "P", "W", "Y", "AA", "R", "T", "U"]


#12/20 
# calstates = ["A", "B", "K"]       ##run0000
# scistates = ["H", "J", "L", "P"]

#calstates = ["B","D", "J"]        ##run0001 
calstates = ["D", "J"]  
scistates = ["F", "N"]
# scistates = ["F"]

#12/21
#calstates = ["A", "B", "I", "O", "AH"]
#scistates = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]



fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)

data = ChannelGroup(fltoday)

defbinsize = 0.5
data.setDefaultBinsize(defbinsize)
data.learnResidualStdDevCut()



ds = data[1]
# ds.plotHist( np.arange(0, 60000, 10), "filtValue", coAddStates=False, states=calstates)   #states=None by default uses all states
# plt.show()
# sfasd
ds.calibrationPlanInit("filtValue")

# #12/15 
# ds.calibrationPlanAddPoint(8950, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10443, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(15546, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(19440, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(25950, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(35590, "FeKAlpha", states=calstates)


# #12/16   
# #run0000              #chan1 
# ds.calibrationPlanAddPoint(9880, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11506, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16957, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20989, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27702, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37443, "FeKAlpha", states=calstates)


# #12/16  
# #run0001          ##chan1 // state F & B 
# ds.calibrationPlanAddPoint(9790, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11430, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16810, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20880, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27500, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37190, "FeKAlpha", states=calstates)



# #12/19        ##run0000
# ds.calibrationPlanAddPoint(8970, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(15377, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(19030, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(25190, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(34610, "FeKAlpha", states=calstates)


#12/20
##run0000       chan 1
# ds.calibrationPlanAddPoint(9793, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11363, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16661, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20668, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27350, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37302, "FeKAlpha", states=calstates)

#run0001
ds.calibrationPlanAddPoint(9740, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(11308, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(16595, "ClKAlpha", states=calstates)
ds.calibrationPlanAddPoint(20580, "KKAlpha", states=calstates)
ds.calibrationPlanAddPoint(27240, "TiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(37180, "FeKAlpha", states=calstates)



#12/21              
# ##run0002             #chan 1
# ds.calibrationPlanAddPoint(8190, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(9550, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(14090, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(17550, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(23270, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(31720, "FeKAlpha", states=calstates)
#ds.plotAvsB("relTimeSec", "filtValue", states=calstates)

#data.learnResidualStdDevCut()
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



external_trigger_filename =  os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_external_trigger.bin")
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)

firstsci = scistates[0]
energies = np.hstack([ds.getAttr("energy", firstsci) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]] for ds in data.values()])

for scistate2 in scistates: 
    energies = np.append(energies,np.hstack([ds.getAttr("energy", scistate2) for ds in data.values()]))
    seconds_after_external_triggers = np.append(seconds_after_external_triggers,np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate2)[0]] for ds in data.values()]))

dat2 = np.vstack((energies,seconds_after_external_triggers))
dat2 = dat2.T 

np.save(datdest+'\\'+str(today)+'_'+str(rn), dat2)

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
#     t_ext = []
#     for i in data:
#         ds = data[i] 
#         energy.extend(list(ds.getAttr('energy', sta)))
#         time.extend(list(ds.getAttr('unixnano', sta)))
#         t_ext.extend(list(ds.seconds_after_external_trigger[ds.getStatesIndicies(states=sta)[0]]))
#     plist = np.array([energy, time, t_ext], dtype=object).T

#     photlist = pd.DataFrame(data=plist, columns=['energy', 'time', 'external_trigger_time'])
#     photlist.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+str('photonlist')+'.csv', index=False)