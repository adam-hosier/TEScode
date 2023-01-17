import mass
import numpy as np
import pylab as plt
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os
import ebit_util
plt.ion()
plt.close("all")

d = "/home/pcuser/data"
date = "20221220"                       #date
rn = "0000"                             #run number, can be found in dastard command 
state1 = "B"
cycle_period = 25 *1e-3
#init_cook = 1
init_cook = 0.980
fl = getOffFileListFromOneFile(os.path.join(d, f"{date}", f"{rn}", 
f"{date}_run{rn}_chan1.off"), maxChans=300)
data = ChannelGroup(fl)             #all of the channels 
data.setDefaultBinsize(0.7)
# 
ds = data[1]                       #selection of individual channel
#ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=state1 )   #states=None by default uses all states



ds.learnDriftCorrection(overwriteRecipe=True, states=state1)
#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
ds.calibrationPlanInit("filtValueDC")
ds.calibrationPlanAddPoint(32261, "FeKAlpha", states=state1)
ds.calibrationPlanAddPoint(17673, "KKAlpha", states=state1)
ds.calibrationPlanAddPoint(23472, "TiKAlpha", states=state1)
ds.calibrationPlanAddPoint(8466, "AlKAlpha", states=state1)
ds.calibrationPlanAddPoint(14265, "ClKAlpha", states=state1)
ds.calibrationPlanAddPoint(9721, "SiKAlpha", states=state1)

#ds.plotAvsB("relTimeSec", "pretriggerMean")
#ds.plotAvsB("relTimeSec", "energyRough")
#plt.ylim(1100,1600)
#ds.learnPhaseCorrection(uncorrectedName="filtValueDC", linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)


#ds.plotHist( np.arange(0,8500,2), "energyRough", coAddStates=False, states=state1)

ds.calibrateFollowingPlan("filtValueDC")
#ds.diagnoseCalibration()
#ds.plotAvsB("filtPhase", "energy")
#ds.plotHist( np.arange(0,14000,1), "energy", coAddStates=False, states=state1)
#ds.plotAvsB("relTimeSec", "energy")
#ds.linefit("FeKAlpha", states=state1)



data.learnDriftCorrection(overwriteRecipe=True, states="B")
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

#states_for_2dhist = 'J'
states_for_2dhist_list = ["H", "J", "L"]
#########

states_for_2dhist = states_for_2dhist_list[0]
energies = np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()])

for states_for_2dhist in states_for_2dhist_list[1:]:
    print(states_for_2dhist)
    energies = np.append(energies,np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()]))
    print([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()])
    seconds_after_external_triggers = np.append(seconds_after_external_triggers,np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()]))

#########
# energies = np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()])
# seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()])

#seconds_after_external_triggers = (seconds_after_external_triggers-init_cook) % cycle_period
#seconds_after_external_triggers = (seconds_after_external_triggers-init_cook)

plt.figure()
A,_,_,_ = plt.hist2d(seconds_after_external_triggers, 
    energies, 
    #bins=(np.arange(0,25/(1e3),1e-4), np.arange(500,2000,3)))
    #bins=(np.arange(0,1000,1), np.arange(500,2000,3)), cmin=1)
    bins=(np.arange(0,1,0.0005), np.arange(500,2000,3)), cmin=1)
print(np.nansum(A))
plt.xlabel("time since external trigger (s)")
plt.ylabel("energy(eV)")
plt.title(f"{data.shortName}, states={states_for_2dhist}")
plt.colorbar()

data2 = np.vstack((energies,seconds_after_external_triggers))
data2 = data2.T
np.save('/home/pcuser/Desktop/TES-GUI/M3 Lifetime/20221220_0000_hjl', data2)
