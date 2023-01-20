import math
import time
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import ebit_util
plt.ion()
d = "C:\\Users\\ahosi\\OneDrive\\Desktop\\tesdata"
date = "20221220"                       #date
rn = "0001"                             #run number, can be found in dastard command 
state1 = "D"                #calibration state
fl = getOffFileListFromOneFile(os.path.join(d, f"{date}", f"{rn}", 
f"{date}_run{rn}_chan1.off"), maxChans=300)
data = ChannelGroup(fl)             #all of the channels 
data.setDefaultBinsize(0.4)

#calibrate and inspect one channel
ds = data[1]                       #selection of individual channel
ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=state1 )   #states=None by default uses all states


ds.learnDriftCorrection(overwriteRecipe=True, states=state1)
#ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
ds.calibrationPlanInit("filtValueDC")
ds.calibrationPlanAddPoint(9740, "AlKAlpha", states=state1)
ds.calibrationPlanAddPoint(11308, "SiKAlpha", states=state1)
ds.calibrationPlanAddPoint(16595, "ClKAlpha", states=state1)
ds.calibrationPlanAddPoint(20580, "KKAlpha", states=state1)
ds.calibrationPlanAddPoint(27240, "TiKAlpha", states=state1)
ds.calibrationPlanAddPoint(37180, "FeKAlpha", states=state1)

ds.plotHist( np.arange(0, 10000, 1), "energyRough", coAddStates=False, states=state1 )   #states=None by default uses all states
#ds.plotAvsB("relTimeSec", "pretriggerMean")
#ds.plotAvsB("relTimeSec", "energyRough")
#plt.ylim(1100,1600)
ds.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)


#ds.plotHist( np.arange(0,8500,2), "energyRough", coAddStates=False, states=state1)

ds.calibrateFollowingPlan("filtValuePC")
ds.diagnoseCalibration()
ds.plotAvsB("filtPhase", "energy")
ds.plotAvsB("pretriggerMean","energy")

data.learnDriftCorrection(overwriteRecipe=True, states=state1)
data.alignToReferenceChannel(ds, "filtValueDC", np.arange(0,40000,30))
data.learnPhaseCorrection(linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
data.calibrateFollowingPlan("filtValuePC")

# def refresh_and_print_if_worked():
#     old = np.sum([len(dslocal) for dslocal in data.values()])
#     data.refreshFromFiles()
#     new = np.sum([len(dslocal) for dslocal in data.values()])
#     if new > old:
#         print(f"refresh worked, now have {new} pulses vs {old} before")
#     else:
#         print(f"refresh FAILED, now have {new} pulses vs {old} before")

#realtime 2d specta, re-run this cell to refresh
external_trigger_filename =  os.path.join(d, f"{date}",f"{rn}", f"{date}_run{rn}_external_trigger.bin")
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)


#states_for_2dhist_list = ["G","K","O","Q", "U", "AD", "AB", "Y","Z", "AH","AO"]
states_for_2dhist_list = ["F"]

states_for_2dhist = states_for_2dhist_list[0]
energies = np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()])
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)[0]] for ds in data.values()])

print(data.shortName)

# for states_for_2dhist in states_for_2dhist_list[1:]:
#     energies = np.append(energies,np.hstack([ds.getAttr("energy", states_for_2dhist) for ds in data.values()]))
#     seconds_after_external_triggers = np.append(seconds_after_external_triggers,np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies(states=states_for_2dhist)] for ds in data.values()]))

data2 = np.vstack((energies,seconds_after_external_triggers))
data2 = data2.T
np.save('C:\\Users\\ahosi\\OneDrive\\Desktop\\testfilem3', data2)