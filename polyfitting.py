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

xdat = np.arange(18, 93)

ydat = [71.3801,
    63.6303,
    57.3997,
    52.2766,
    47.9876,
    44.3421,
    41.2039,
    38.4733,
    36.0745,
    33.9499,
    32.0547,
    30.3529,
    28.8161,
    27.4209,
    26.1485,
    24.9827,
    23.9108,
    22.9213,
    22.0052,
    21.1541,
    20.3613,
    19.6209,
    18.9276,
    18.2772,
    17.6655,
    17.0889,
    16.5446,
    16.0296,
    15.5417,
    15.0786,
    14.6387,
    14.2197,
    13.8205,
    13.4396,
    13.0758,
    12.7272,
    12.3938,
    12.0739,
    11.767,
    11.4721,
    11.1886,
    10.9156,
    10.6527,
    10.3994,
    10.1553,
    9.9191,
    9.6911,
    9.4704,
    9.2566,
    9.051,
    8.8496,
    8.6562,
    8.4693,
    8.2856,
    8.1085,
    7.9362,
    7.7696,
    7.6065,
    7.4484,
    7.2948,
    7.145,
    6.9994,
    6.8579,
    6.7198,
    6.5854,
    6.4547,
    6.3269,
    6.2033,
    6.0837,
    5.9656,
    5.8512,
    5.7393,
    5.6358,
    5.5251,
    5.4296]



def fun1(x, a, b, c, d, e, f, g):
    return a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2 +f*x + g

def fun2(x, a, b, c, d):
    return b + a/x + c/(x**2) + d/ (x**3)


model2 = Model(fun2)
params = Parameters() 
params.add('a', value = 0, vary=True)
params.add('b', value = 0, vary=True)
params.add('c', value = 0, vary=True)
params.add('d', value = 0, vary=True)
# params.add('e', value = 1)
# params.add('f', value = 1)
# params.add('g', value = 1)
fit2 = model2.fit(ydat, params, x=xdat, weights=None)
params.update(fit2.params)
xplot = np.linspace(np.min(xdat), np.max(xdat), num=1000)
yplot = model2.eval(params=params, x=xplot)


params.pretty_print()
plt.figure() 
plt.plot(xplot, yplot, label='fit', c='r')
plt.scatter(xdat, ydat)
plt.legend()
plt.ylabel('wavelength [nm]')
plt.xlabel('Z')
plt.show()
plt.close() 