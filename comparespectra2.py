import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util

d = "/home/pcuser/data"
today = "20221221"
todayrn = "0002"

yesterday = "20221219"
yesterdayrn = "0000"

todaycal = "AH"
yesterdaycal = "C"

fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{todayrn}", 
f"{today}_run{todayrn}_chan1.off"), maxChans=300)

flyes = getOffFileListFromOneFile(os.path.join(d, f"{yesterday}", f"{yesterdayrn}", 
f"{yesterday}_run{yesterdayrn}_chan1.off"), maxChans=300)


datatod = ChannelGroup(fltoday)
datayes = ChannelGroup(flyes)

datatod.setDefaultBinsize(0.5)
datayes.setDefaultBinsize(0.7)

# ##### Yesterday 
# ds = datayes[1]
# #ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=yesterdaycal )   #states=None by default uses all states

# ds.learnDriftCorrection(overwriteRecipe=True, states=yesterdaycal)
# #ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
# ds.calibrationPlanInit("filtValueDC")

# #ds.calibrationPlanAddPoint(34594, "FeKAlpha", states=yesterdaycal)
# #ds.calibrationPlanAddPoint(19010, "KKAlpha", states=yesterdaycal)
# #ds.calibrationPlanAddPoint(25180, "TiKAlpha", states=yesterdaycal)
# ds.calibrationPlanAddPoint(8980, "AlKAlpha", states=yesterdaycal)
# ds.calibrationPlanAddPoint(15387, "ClKAlpha", states=yesterdaycal)
# ds.calibrationPlanAddPoint(10450, "SiKAlpha", states=yesterdaycal)

# ds.calibrateFollowingPlan("filtValueDC")

# datayes.learnDriftCorrection(overwriteRecipe=True, states=yesterdaycal)
# datayes.alignToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
# datayes.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
# datayes.calibrateFollowingPlan("filtValuePC")




#### Today 
ds2 = datatod[1] 

ds2.learnDriftCorrection(overwriteRecipe=True, states=todaycal)
ds2.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states="AH")   #states=None by default uses all states
ds2.calibrationPlanInit("filtValueDC")
ds2.calibrationPlanAddPoint(31630, "FeKAlpha", states=todaycal)
ds2.calibrationPlanAddPoint(17550, "KKAlpha", states=todaycal)
ds2.calibrationPlanAddPoint(23260, "TiKAlpha", states=todaycal)
ds2.calibrationPlanAddPoint(8170, "AlKAlpha", states=todaycal)
ds2.calibrationPlanAddPoint(14030, "ClKAlpha", states=todaycal)
ds2.calibrationPlanAddPoint(9540, "SiKAlpha", states=todaycal)

ds2.calibrateFollowingPlan("filtValueDC")



datatod.learnDriftCorrection(overwriteRecipe=True, states=todaycal)
datatod.alignToReferenceChannel(ds2, "filtValueDC", np.arange(0,40000,30))
datatod.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
datatod.calibrateFollowingPlan("filtValuePC")


Ndstates = ["H", "W", "Y", "AA"]

# datayes.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=Ndstates)

datatod.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=["E","G","K","M","R"])


