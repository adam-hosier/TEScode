import mass
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os
import h5py
import pylab as plt
import numpy as np
import lmfit
plt.ion()

pulse_files = "/home/pcuser/data/20221214/0001/20221214_run0001_chan*.ljh"
noise_files = "/home/pcuser/data/20221214/0000/20221214_run0000_chan*.ljh"



data_plain = mass.TESGroup(filenames = pulse_files,
        noise_filenames = noise_files,
        overwrite_hdf5_file=True)
data_plain.summarize_data(pretrigger_ignore_microsec=500)
data_plain.auto_cuts()
data_plain.compute_noise(max_excursion=300)
data_plain.compute_ats_filter(f_3db=10e3)
data_plain.filter_data()

ds = data_plain.channel[5] 
for ds in data_plain:
        # ds.calibrate("p_filt_value",
        # line_names=["KKAlpha", "ClKAlpha", "FeKAlpha"])
        dso = ds.toOffStyle()
        dso.plotHist(np.arange(100,43000,40), "p_filt_value")