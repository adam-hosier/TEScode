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

todaycal = ["A", "B", "I", "O", "AH"]


datdest = 'C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022'

fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan1.off"), maxChans=300)


data = ChannelGroup(fltoday)

defbinsize = 0.1
data.setDefaultBinsize(defbinsize)



ds = data[1]
ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=todaycal)   #states=None by default uses all states
# plt.show()
# asdfas


#ds.learnTimeDriftCorrection(states=todaycal)
ds.learnDriftCorrection(overwriteRecipe=True, states=todaycal)

#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
ds.calibrationPlanInit("filtValueDC")

#12/15 
# ds.calibrationPlanAddPoint(8950, "AlKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(15560, "ClKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(19440, "KKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(25950, "TiKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(35570, "FeKAlpha", states=todaycal)


# #12/16 
# ds.calibrationPlanAddPoint(8660, "AlKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(10049, "SiKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(14785, "ClKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(18292, "KKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(24085, "TiKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(32678, "FeKAlpha", states=todaycal)


# #12/19 
# ds.calibrationPlanAddPoint(8970, "AlKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(10460, "SiKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(15377, "ClKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(19030, "KKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(25190, "TiKAlpha", states=todaycal)
# ds.calibrationPlanAddPoint(34610, "FeKAlpha", states=todaycal)


#12/21
ds.calibrationPlanAddPoint(8190, "AlKAlpha", states=todaycal)
ds.calibrationPlanAddPoint(9550, "SiKAlpha", states=todaycal)
ds.calibrationPlanAddPoint(14090, "ClKAlpha", states=todaycal)
ds.calibrationPlanAddPoint(17550, "KKAlpha", states=todaycal)
ds.calibrationPlanAddPoint(23270, "TiKAlpha", states=todaycal)
ds.calibrationPlanAddPoint(31720, "FeKAlpha", states=todaycal)


#ds.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
# ds.calibrateFollowingPlan("filtValuePC")

ds.calibrateFollowingPlan("filtValueDC")

data.learnDriftCorrection(overwriteRecipe=True, states=todaycal)

data.alignToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
#data.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
data.calibrateFollowingPlan("filtValueDC")


scistates = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]

#data.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=scistates)

for sta in scistates:
    histall = np.array(data.hist(np.arange(0, 14000, defbinsize), "energy", states=sta))
    stadat = pd.DataFrame(data=histall)
    #stadat = stadat.transpose
    stadat.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+'.csv', index=False)

    energy = []
    time = []
    for i in data:
        ds = data[i] 
        energy.extend(list(ds.getAttr('energy', sta)))
        time.extend(list(ds.getAttr('unixnano', sta)))
    plist = np.array([energy, time], dtype=object).T

    photlist = pd.DataFrame(data=plist, columns=['energy', 'time'])
    photlist.to_csv(datdest + '/' + str(today) + '_' + str(rn) + '_' + str(sta)+str('photonlist')+'.csv', index=False)

#ljh2off C:\Users\ahosi\OneDrive\Desktop\tesdata\20221221\0002\20221221_run0002_chan1.ljh C:\Users\ahosi\OneDrive\Desktop\tesdata\20221214\0001\20221214_run0001_model.hdf5 C:\Users\ahosi\OneDrive\Desktop\tesdata\20221221\0002_newoff -r

#plt.show()

