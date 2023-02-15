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

t = [0,
    0.8334,
    4,
    17.5,
    22.75,
    23.5]

p = [94,
    79,
    62,
    38,
    37,
    36]

def mod1(x, A, tau, B):
    return A*np.exp(-x/tau)+B


model1 = Model(mod1)
params = Parameters() 
params.add('A', value=40)
params.add('tau', value=2)
params.add('B',value=60)

fit = model1.fit(p, params, x=t, weights=None)
params.update(fit.params)

xplot = np.linspace(np.min(t), np.max(t), num=1000)
yplot = model1.eval(params=params, x=xplot)

plt.figure() 
plt.plot(xplot, yplot, label='fit', c='r')
plt.scatter(t, p, label='pressure data mTorr')
plt.legend()
plt.ylabel('pressure (mTorr)')
plt.xlabel('time [hrs]')
plt.show()
plt.close() 

