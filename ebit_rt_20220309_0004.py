import mass
import numpy as np
import pylab as plt
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os
import ebit_util
plt.ion()
plt.close("all")

d = "/home/pcuser/data"
date = "20220309"                       #date
rn = "0004"                             #run number, can be found in dastard command 
state1 = "D"
fl = getOffFileListFromOneFile(os.path.join(d, f"{date}", f"{rn}", 
f"{date}_run{rn}_chan1.off"), maxChans=300)
data = ChannelGroup(fl)             #all of the channels 
data.setDefaultBinsize(0.7)
# 
ds = data[17]                       #selection of individual channel
ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=state1 )   #states=None by default uses all states


ds.learnDriftCorrection(overwriteRecipe=True, states=state1)
#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
ds.calibrationPlanInit("filtValueDC")
#ds.calibrationPlanAddPoint(37150, "CuKAlpha", states=state1)
ds.calibrationPlanAddPoint(34225, "FeKAlpha", states=state1)
ds.calibrationPlanAddPoint(37150, "FeKBeta", states=state1)
#ds.calibrationPlanAddPoint(2661, "TiKAlpha", states=state1)
ds.calibrationPlanAddPoint(9048, "AlKAlpha", states=state1)
ds.calibrationPlanAddPoint(7690, "MgKAlpha", states=state1)
#ds.calibrationPlanAddPoint(11637, "SiKAlpha", states=state1)

#ds.plotAvsB("relTimeSec", "pretriggerMean")
#ds.plotAvsB("relTimeSec", "energyRough")
#plt.ylim(1100,1600)
ds.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)


#ds.plotHist( np.arange(0,8500,2), "energyRough", coAddStates=False, states=state1)

ds.calibrateFollowingPlan("filtValuePC")
ds.diagnoseCalibration()
#ds.plotAvsB("filtPhase", "energy")
#ds.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=state1)
#ds.plotAvsB("relTimeSec", "energy")
#ds.linefit("FeKAlpha", states=state1)



data.learnDriftCorrection(overwriteRecipe=True, states="D")
data.alignToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
data.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
data.calibrateFollowingPlan("filtValuePC")
#data.linefit("FeKAlpha", states="F")
#data.qualityCheckLinefit("AlKAlpha", worstAllowedFWHM=5, states="F")
#data.linefit("AlKAlpha", states="F")
#data.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states="F")



external_trigger_filename =  os.path.join(d, f"{date}",f"{rn}", f"{date}_run{rn}_external_trigger.bin")
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)

states_for_2dhist = "H"
energies = np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()])

plt.figure()
plt.hist2d(seconds_after_external_triggers, 
    energies, 
    bins=(np.arange(0,3,0.01), np.arange(500,5000,3)))
plt.xlabel("time since external trigger (s)")
plt.ylabel("energy(eV)")
plt.title(f"{data.shortName}, states={states_for_2dhist}")
plt.colorbar()

