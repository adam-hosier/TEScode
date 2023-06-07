import mass 
import numpy as np 
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
#import ebit_util
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import binned_statistic_2d
from scipy import stats
from fit_utils import MultiPeakGaussian
from lmfit import minimize, Parameters, report_fit, Model, Parameter
from lmfit.models import GaussianModel
from lmfit.models import SplitLorentzianModel 


n = [4, 5, 8, 9]
E = [1203.89,
    1545.7,
    1899,
    1948]

#2134 IE  Ni-like Nd    
#2224 IE Co-like Nd

### Ni-like series
E = [930.11,
    588.3,
    235,
    186]

# n = [4, 5, 6]
# E = [930.11,
#     588.3,
#     422]

ntest = [4, 5, 6, 7, 8, 9, 10]

# n = [4, 5, 6]
# E2 = [983, 618, 491]

def mod2(x, A, B, C):
    return C + A/(x+B)**2



model2 = Model(mod2)
params = Parameters() 
params.add('A', value = 14816)
params.add('B', value = 0.001)
params.add('C', value = 0, vary=False)
fit2 = model2.fit(E, params, x=n, weights=None)
params.update(fit2.params)
xplot = np.linspace(np.min(n), np.max(n), num=1000)
yplot = model2.eval(params=params, x=xplot)


params.pretty_print()
plt.figure() 
plt.plot(xplot, yplot, label='fit', c='r')
plt.scatter(n, E, label='Energy')
plt.legend()
plt.ylabel('ionization energy - photon energy [eV]')
plt.xlabel('n (principal quantum number)')
plt.show()
plt.close() 