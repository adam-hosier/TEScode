import mass
import numpy as np
import pylab as plt

def get_external_triggers(filename, good_only):
    f = open(filename,"rb")
    f.readline() # discard comment line
    external_trigger_rowcount = np.fromfile(f,"int64")
    if good_only:
        external_trigger_rowcount = external_trigger_rowcount[get_good_trig_inds(external_trigger_rowcount)]
    return external_trigger_rowcount

def get_good_trig_inds(external_trigger_rowcount, plot=False):
    # ignore triggers too close to each other
    d = np.diff(external_trigger_rowcount)
    median_diff = np.median(np.diff(external_trigger_rowcount))
    good_inds = np.where(d > median_diff/2)[0]
    if plot:
        plt.figure()
        plt.plot(d,".",label="all")
        plt.plot(good_inds, d[good_inds],".",label="good")
        plt.legend()
    return good_inds

def calc_external_trigger_timing(self, external_trigger_rowcount):
    nRows = self.offFile.header["ReadoutInfo"]["NumberOfRows"]
    rowcount = self.offFile["framecount"] * nRows
    rows_after_last_external_trigger, rows_until_next_external_trigger = \
        mass.core.analysis_algorithms.nearest_arrivals(rowcount, external_trigger_rowcount)
    self.rowPeriodSeconds = self.offFile.framePeriodSeconds/float(nRows)
    self.rows_after_last_external_trigger = rows_after_last_external_trigger
    self.rows_until_next_external_trigger = rows_until_next_external_trigger
    self.seconds_after_external_trigger = rows_after_last_external_trigger*self.rowPeriodSeconds
    self.seconds_until_next_external_trigger = rows_until_next_external_trigger*self.rowPeriodSeconds
