import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
import ebit_util

d = os.path.dirname(os.path.abspath(__file__))
today = "20230630"
rn = "0001"
fltoday = getOffFileListFromOneFile(os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_chan4.off"), maxChans=300)
external_trigger_filename =  os.path.join(d, f"{today}", f"{rn}", 
f"{today}_run{rn}_external_trigger.bin")
data = ChannelGroup(fltoday)
data.setDefaultBinsize(0.5)
ds = data[1]
external_trigger_rowcount = ebit_util.get_external_triggers(external_trigger_filename, good_only=True)
for ds in data.values():
    ebit_util.calc_external_trigger_timing(ds, external_trigger_rowcount)
seconds_after_external_triggers = np.hstack([ds.seconds_after_external_trigger[ds.getStatesIndicies()[0]] for ds in data.values()])
seconds_until_next_external_trigger = np.hstack([ds.seconds_until_next_external_trigger[ds.getStatesIndicies()[0]] for ds in data.values()])
rows_after_last_external_trigger = np.hstack([ds.rows_after_last_external_trigger[ds.getStatesIndicies()[0]] for ds in data.values()])
rows_until_next_external_trigger = np.hstack([ds.rows_until_next_external_trigger[ds.getStatesIndicies()[0]] for ds in data.values()])
nanotime = np.hstack([ds.unixnano[ds.getStatesIndicies()[0]] for ds in data.values()])

timediff = []

for i in range(len(nanotime)-1):
    timediff.append(nanotime[i+1]-nanotime[i])

rowdiff = []
for i in range(len(rows_after_last_external_trigger)-1):
    rowdiff.append(rows_after_last_external_trigger[i+1]-rows_after_last_external_trigger[i])

plt.hist(timediff, bins=100)
plt.title('time stamp difference histogram')
plt.show()

plt.hist(rowdiff, bins=100)
plt.title('row count since last external trigger difference histogram')
plt.show()