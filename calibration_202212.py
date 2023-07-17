import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util
import pandas as pd
#plt.ion()
d = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata'
#d = 'C:\\data\\tesdata'
today = "20221221"
rn = "0002"
# datdest = 'C:\\data\\TES_Spectra_1eVbin'
# datdest = 'C:\\data\\TES_newcal'
datdest = 'C:\\data\\TimingData'

savedat = False

#12/14 
#calstates = ["B", "C"]

#12/15 ##run0001        
# calstates = ["A", "B", "I", "M", "S", "W", "AF"]
# scistates = ["G", "K", "O", "Q", "U" ,"Y", "Z", "AB", "AD", "AH", "AO"]
#scistates = calstates


#12/16
# calstates = ["A", "B", "F", "J", "N", "R"]    ##run0000
# #scistates = ["E", "H", "L", "P"]
# scistates = calstates

# calstates = ["B", "F"]   ##run0001
# #scistates = ["D"]
# scistates = calstates


#12/19          ## run 0000
# #calstates = ["A", "B", "C", "D", "AC"]
# #calstates = ["B", "C", "D", "AC"]
# calstates = ["AC"]
# scistates = ["H", "K", "P", "W", "Y", "AA", "R", "T", "U"]
# #cistates = ["H", "K", "W", "Y", "AA", "R", "T", "U"]
# # scistates = calstates


#12/20 
# calstates = ["A", "B", "K"]       ##run0000
# #scistates = ["H", "J", "L", "P"]
# scistates = calstates


# calstates = ["B","D", "J"]        ##run0001 
# # #calstates = ["D", "J"]  
# scistates = ["F", "N"]
# #scistates = ["F"]
# scistates = calstates 


#12/21      ##run0002
calstates = ["A", "I", "O", "AH"]
#calstates = ["O"]
scistates = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]
#scistates = calstates 


fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)

data = ChannelGroup(fltoday)

defbinsize = 0.5
data.setDefaultBinsize(defbinsize)
#data.learnResidualStdDevCut()


plotbinSize = 1
ds = data[1]

#ds.plotHist( np.arange(0, 60000, 10), "filtValue", coAddStates=False, states=calstates)   #states=None by default uses all states


#         
# for k in range(1, len(data)):
#     try:
#         ds = data[k]
#         #ds.plotHist( np.arange(0, 60000, 10), "filtValue", coAddStates=False, states='B')   #states=None by default uses all states
#         ds.plotHist( np.arange(0, 60000, 10), "filtValue", coAddStates=False, states=calstates)
#         plt.show() 
#         plt.close() 
#     except:
#         pass 

ds.calibrationPlanInit("filtValue")

# #12/15 
# ds.calibrationPlanAddPoint(8950, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(10443, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(15546, "ClKAlpha", states=calstates)
# #ds.calibrationPlanAddPoint(16660, "ClKBeta", states=calstates)
# ds.calibrationPlanAddPoint(19440, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(25950, "TiKAlpha", states=calstates)
# # ds.calibrationPlanAddPoint(28175, "TiKBeta", states=calstates)
# ds.calibrationPlanAddPoint(35590, "FeKAlpha", states=calstates)
# # ds.calibrationPlanAddPoint(38770, "FeKBeta", states=calstates)


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
# #ds.calibrationPlanAddPoint(16436, "ClKBeta", states=calstates)
# ds.calibrationPlanAddPoint(19030, "KKAlpha", states=calstates)
# #ds.calibrationPlanAddPoint(25190, "TiKAlpha", states=calstates)
# #ds.calibrationPlanAddPoint(27375, "TiKBeta", states=calstates)
# #ds.calibrationPlanAddPoint(34610, "FeKAlpha", states=calstates)
# #ds.calibrationPlanAddPoint(37719, "FeKBeta", states=calstates)

#12/20
# ##run0000       chan 1
# ds.calibrationPlanAddPoint(9793, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11363, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16661, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20668, "KKAlpha", states=calstates)
# #ds.calibrationPlanAddPoint(27350, "TiKAlpha", states=calstates)
# #ds.calibrationPlanAddPoint(37302, "FeKAlpha", states=calstates)

#run0001
# ds.calibrationPlanAddPoint(9740, "AlKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(11308, "SiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(16595, "ClKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(20580, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(27240, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(37180, "FeKAlpha", states=calstates)



#12/21              
# ##run0002             #chan 1
ds.calibrationPlanAddPoint(8190, "AlKAlpha", states=calstates)
ds.calibrationPlanAddPoint(9550, "SiKAlpha", states=calstates)
ds.calibrationPlanAddPoint(14070, "ClKAlpha", states=calstates)
#ds.calibrationPlanAddPoint(15080, "ClKBeta", states=calstates)
ds.calibrationPlanAddPoint(17570, "KKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(23279, "TiKAlpha", states=calstates)
# ds.calibrationPlanAddPoint(31720, "FeKAlpha", states=calstates)



# ds.plotAvsB("relTimeSec", "filtValue", states=calstates)
# plt.grid()
# plt.show() 
# asdf

data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
#data.alignToReferenceChannel(referenceChannel=ds, binEdges=np.arange(0,60000,10), attr="filtValue", states=calstates)
# Phase, drift, and time drift correct on pulses in the rough 400-2500 eV range
#data.cutAdd("Cut1", lambda energyRough: np.logical_and(energyRough > 400, energyRough < 3000))
data.learnPhaseCorrection(uncorrectedName="filtValue", correctedName = "filtValuePC", states=calstates)

data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValuePC", correctedName = "filtValuePCDC", states=calstates)
#data.learnDriftCorrection(indicatorName="pretriggerMean", uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC",)
#data.learnDriftCorrection(uncorrectedName="filtValue", correctedName = "filtValueDC", states=calstates, cutRecipeName="cutForLearnDC",)

data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates) 
#data.learnTimeDriftCorrection(indicatorName="relTimeSec", uncorrectedName="filtValuePCDC", correctedName = "filtValuePCDCTC", states=calstates, _rethrow=True, cutRecipeName="cutForLearnDC",) 

data.calibrateFollowingPlan("filtValuePCDCTC", calibratedName="energy", overwriteRecipe=True)      
#data.calibrateFollowingPlan("filtValuePCDC", calibratedName="energy", overwriteRecipe=True)      

#ds.diagnoseCalibration() 


# ddest = "C:\\Users\\ahosi\OneDrive\Desktop\\TES_Calibration_Lines\\CalibrationXY"

# llist = ['AlKAlpha', 'SiKAlpha', 'ClKAlpha', 'KKAlpha', 'TiKAlpha', 'FeKAlpha']
# for li in llist:
#     temp = data.linefit(str(li), states=calstates)

#     xdat = temp.userkws['bin_centers']
#     ydat = temp.data 
#     aldat = np.array((xdat, ydat))
#     aldat = aldat.T 
#     ndat = pd.DataFrame(data=aldat, columns=['Energy [eV]', 'Counts'])
#     ndat.to_csv(ddest + '/' + str(today) + '_' + str(rn) + '_' + str(li)+'_'+str('sum')+'.csv', index=False)

plt.close() 
#ddest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\TES_Calibration_Lines\\20221221'

#ddest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\TES_Calibration_Lines\\SDr_Test'

# ddest = datdest
# allstates = calstates + scistates
# newhist = data.plotHist( np.arange(200,10000,1), "energy", coAddStates=False, states=allstates)

# linep = plt.gca() 
# conames = []

# for k in range(0, len(linep.lines)):
#     #temp = data.linefit(str(li), states=calstates)
#     temp = linep.lines[k]
#     if k==0: 
#         aldat = np.array(temp.get_xydata())
#         conames.append(str(today)+'_'+str(rn)+'_'+str(temp.get_label())+'_Energy')
#     else: 
#         len1 = len(temp.get_xydata())
#         aldat = np.hstack((aldat, temp.get_xydata()[:,1].reshape((len1),1)))

#     calna = ''.join(calstates)
#     #conames.append(str(today)+'_'+str(rn)+'_'+str(temp.get_label())+'_Energy')
#     conames.append(str(today)+'_'+str(rn)+'_'+str(calna)+'cal'+'_'+str(temp.get_label())+'_Counts')
#     #aldat = aldat.T 
#     # ndat = pd.DataFrame(data=aldat, columns=['Energy [eV]', 'Counts'])
#     # ndat.to_csv(ddest + '/' + str(today) + '_' + str(rn) + '_' + str(temp.get_label())+'_'+str('testDR')+'.csv', index=False)

# ndat = pd.DataFrame(data=aldat, columns=conames)
# #ndat.to_csv(ddest + '/' + str('AlSiClKTiFe_only')+'.csv', index=False)

# newname = ''.join(calstates)
# ndat.to_csv(ddest + '/' +str(today)+'_'+str(rn)+'_'+newname+ str('_cal')+'.csv', index=False)


# for l in llist: 
#     temp = data.linefit(str(l), states=calstates)
#     print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
#     print(str(l))

#     for k in temp.params.keys():
#         print(str(k),' :',temp.params[str(k)].value ,' +/- ' ,temp.params[str(k)].stderr)


external_trigger_filename =  os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_external_trigger.bin")
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)



firstsci = scistates[0]
# #firstsci = "F"
energies = np.hstack([ds.getAttr("energy", firstsci) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]] for ds in data.values()])

#np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=firstsci)[0]][ds.getAttr("cutResidualStdDev", firstsci)] for ds in data.values()])

for scistate in scistates: 
    # energies = np.append(energies,np.hstack([ds.getAttr("energy", scistate) for ds in data.values()]))
    # seconds_after_external_triggers = np.append(seconds_after_external_triggers,np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate)[0]] for ds in data.values()]))

    energies = np.hstack([ds.getAttr("energy", scistate) for ds in data.values()])
    seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=scistate)[0]] for ds in data.values()])


    dat2 = np.vstack((energies,seconds_after_external_triggers))
    dat2 = dat2.T 
    #if savedat: 
    np.save(datdest+'\\'+str(today)+'_'+str(rn)+'_'+str(scistate), dat2)

#data.plotHist(np.arange(0,10000,1),"energy", states=scistates, coAddStates=False)


# data.plotHist(np.arange(0,60000,10),"energy", states=scistates, coAddStates=True)
# ds.plotHist(np.arange(0,60000,10),"filtValue", states=scistates, coAddStates=True)
#data.plotHist( np.arange(0,10000,1), "energy", coAddStates=True, states=calstates)
# plt.show()
# ds.diagnoseCalibration()


# for sta in scistates:
#     histall = np.array(data.hist(np.arange(200, 9000, plotbinSize), "energy", states=sta))
#     histall = histall.T
#     stadat = pd.DataFrame(data=histall)
#     #stadat = stadat.transpose
    
#     stadat.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+'.csv', index=False)

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
    
    # photlist.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+str('photonlist')+'.csv', index=False)