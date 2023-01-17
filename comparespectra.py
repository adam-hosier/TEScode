import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util

d = "/home/pcuser/data"
today = "20221219"
#todayrn = "0000"
todayrn = "0000"

yesterday = "20221219"
yesterdayrn = "0000"

todaycal = "C"
yesterdaycal = "C"

fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{todayrn}", 
f"{today}_run{todayrn}_chan1.off"), maxChans=300)

flyes = getOffFileListFromOneFile(os.path.join(d, f"{yesterday}", f"{yesterdayrn}", 
f"{yesterday}_run{yesterdayrn}_chan1.off"), maxChans=300)


datatod = ChannelGroup(fltoday)
datayes = ChannelGroup(flyes)

datatod.setDefaultBinsize(0.7)
datayes.setDefaultBinsize(1)

##### Yesterday 
ds = datayes[1]
ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=yesterdaycal)   #states=None by default uses all states

ds.learnDriftCorrection(overwriteRecipe=True, states=yesterdaycal)
#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
ds.calibrationPlanInit("filtValueDC")

ds.calibrationPlanAddPoint(34594, "FeKAlpha", states=yesterdaycal)
ds.calibrationPlanAddPoint(19010, "KKAlpha", states=yesterdaycal)
ds.calibrationPlanAddPoint(25180, "TiKAlpha", states=yesterdaycal)
ds.calibrationPlanAddPoint(8980, "AlKAlpha", states=yesterdaycal)
ds.calibrationPlanAddPoint(15387, "ClKAlpha", states=yesterdaycal)
ds.calibrationPlanAddPoint(10450, "SiKAlpha", states=yesterdaycal)

ds.calibrateFollowingPlan("filtValueDC")

datayes.learnDriftCorrection(overwriteRecipe=True, states=yesterdaycal)
datayes.alignToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
#datayes.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
datayes.calibrateFollowingPlan("filtValueDC")

Ndstates = ["AA","W","Y","M","P"]

datayes.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=Ndstates)
#plt.show()


# #### Today 
# ds2 = datatod[1] 

# ds2.learnDriftCorrection(overwriteRecipe=True, states=todaycal)
# #ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
# ds2.calibrationPlanInit("filtValueDC")
# # ds2.calibrationPlanAddPoint(32261, "FeKAlpha", states=todaycal)
# # ds2.calibrationPlanAddPoint(17673, "KKAlpha", states=todaycal)
# # ds2.calibrationPlanAddPoint(23472, "TiKAlpha", states=todaycal)
# ds2.calibrationPlanAddPoint(8466, "AlKAlpha", states=todaycal)
# ds2.calibrationPlanAddPoint(14265, "ClKAlpha", states=todaycal)
# ds2.calibrationPlanAddPoint(9721, "SiKAlpha", states=todaycal)

# ds2.calibrateFollowingPlan("filtValueDC")



# datatod.learnDriftCorrection(overwriteRecipe=True, states=todaycal)
# datatod.alignToReferenceChannel(ds2, "filtValueDC", np.arange(0,40000,30))
# datatod.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
# datatod.calibrateFollowingPlan("filtValuePC")


# Ndstates = ["AA", "W", "Y", "M", "P"]

# datayes.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=Ndstates)
# datayes.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=["AA", "W", "Y", "M", "P"])

