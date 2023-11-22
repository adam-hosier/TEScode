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

##################
# Co-like series 
## f shell HIGHER E series 3d9 nf1
#2224 IE Co-like Nd

n = [4,5,6]
TE = [1172.29,1522.32,1711.60]

# ## f shell HIGHER E series 3d9 nf1
# n = [4,5,6]
# TE = [1203.89, 1546, 1733.94]



# n = [4, 5, 8, 9]
# E = [1203.89,
#     1545.7,
#     1899,
#     1948]

#2134 IE  Ni-like Nd    

######################
# Ni-like series
# #f shell lower energy series 3d8 nf1
# n = [4,5]
# TE = [1209,1580.66]
# # E = [930.11,
# #     588.3,
# #     235,
# #     186]

# Ni-like series
#f shell lower energy series 3d8 nf1
n = [4,5]
TE = [1241,1606]

E = [2134-i for i in TE]

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
xplot = np.linspace(np.min(n), np.max(ntest), num=1000)
yplot = model2.eval(params=params, x=xplot)

yeval = model2.eval(params=params, x=ntest)
yevaln = [2134-i for i in yeval]
print(yevaln)
params.pretty_print()
plt.figure() 
plt.plot(xplot, yplot, label='Fit', c='r')
plt.scatter(n, E)
plt.scatter(ntest, yeval, label='Extracted Energies', c='g')
plt.legend()
plt.ylabel('Ionization energy - Photon energy [eV]')
plt.xlabel('n (Principal Quantum Number)')
plt.show()
plt.close() 
