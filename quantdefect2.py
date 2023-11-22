import numpy as np 
import pylab as plt 
import os
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import binned_statistic_2d
from scipy import stats
from lmfit import minimize, Parameters, report_fit, Model, Parameter
from lmfit.models import GaussianModel
from lmfit.models import SplitLorentzianModel 
#2134 IE  Ni-like Nd  (33)
#2224 IE Co-like Nd (34)
IE = 2224
##################
# Ni-like series 
## f shell lower E series 3d9 nf1


# n = [4,5,6]
# TE = [1158,1522.32,1711.60]



# ## f shell HIGHER E series 3d9 nf1
# n = [4,5, 6, 7, 8]
# TE = [1203.89, 1546, 1733.94, 1847.57, 1920.77]
# n = [5, 6, 7, 8]
# TE = [ 1546, 1733.94, 1847.57, 1920.77]



# n = [4, 5, 8, 9]
# E = [1203.89,
#     1545.7,
#     1899,
#     1948]

  

######################
# Co-like series
# #f shell lower energy series 3d8 nf1
# n = [4,5,6]
# TE = [1241,1580.66,1774.46]
n = [5,6]
TE = [1580.66,1774.46]
n = [4,5]
TE = [1241,1580.66]

# Co-like series
#f shell higher energy series 3d8 nf1
# n = [4,5]
# TE = [1241,1606]


#p shell lower energy 
# n = [4,5]
# TE = [931.36, 1416.49]

E = [IE-i for i in TE]
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

yeval = model2.eval(params=params, x=ntest)
yevaln = [IE-i for i in yeval]
print(yevaln)
params.pretty_print()
plt.figure() 
plt.plot(xplot, yplot, label='fit', c='r')
plt.scatter(n, E, label='Energy')
plt.legend()
plt.ylabel('ionization energy - photon energy [eV]')
plt.xlabel('n (principal quantum number)')
plt.show()
plt.close() 
