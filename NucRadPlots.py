import numpy as np 
import pylab as plt 
import os 
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import binned_statistic_2d
from scipy import stats
#from fit_utils import MultiPeakGaussian
from lmfit import minimize, Parameters, report_fit, Model, Parameter
from lmfit.models import GaussianModel
from lmfit.models import SplitLorentzianModel 
from matplotlib.ticker import MaxNLocator

from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
W1999A = [182,
    183,
    184,
    186]

W1999R = [5.355,
    5.33,
    5.364,
    5.3808]

W1999dR = [0.002,
    0.15,
    0.002,
    0.0062]

Os1999A = [186,
    187,
    188,
    190,
    192]

Os1999R = [5.387,
    5.4001,
    5.36,
    5.4062,
    5.41]

Os1999dR = [0.0072,
    0.0013,
    0.15,
    0.0014,
    0.002]

Ir1999A = [191, 193]

Ir1999R = [5.39,
    5.41]

Ir1999dR = [0.15,
    0.15]

Pt1999A = [194,
    195,
    196,
    198]

Pt1999R = [5.4043,
    5.427,
    5.3801,
    5.441]

Pt1999dR = [0.0195,
    0.0069,
    0.0115,
    0.0063]

######################################

W2004A = [180,
    182,
    183,
    184,
    186]

W2004R = [5.3493,
5.3566,
5.33,
5.367,
5.3759]

W2004dR = [0.0023,
0.0017,
0.15,
0.0017,
0.0019]

Re2004A = [185,
187]

Re2004R = [5.3286,
5.3386
]

Re2004dR = [0.0126,
0.0127
]

Os2004A = [184,
186,
187,
188,
189,
190,
192
]

Os2004R = [5.382,
5.3908,
5.3934,
5.3994,
5.36,
5.4061,
5.4126
]

Os2004dR = [0.0023,
0.0016,
0.0017,
0.0011,
0.15,
0.0009,
0.0011
]

Ir2004A = [182,
183,
184,
185,
186,
187,
188,
189,
191,
193
]

Ir2004R = [5.3809,
5.3857,
5.3881,
5.3899,
5.3943,
5.3881,
5.389,
5.3938,
5.3981,
5.4019
]

Ir2004dR = [0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
]

Pt2004A = [183,
184,
185,
186,
187,
188,
189,
190,
191,
192,
193,
194,
195,
196,
198
]

Pt2004R = [5.4028,
5.4155,
5.4046,
5.4074,
5.4059,
5.4068,
5.4068,
5.4117,
5.411,
5.4181,
5.4201,
5.4247,
5.4278,
5.4315,
5.4403
]

Pt2004dR = [0.0037,
0.0029,
0.0039,
0.0041,
0.0035,
0.0035,
0.0037,
0.0032,
0.0032,
0.0029,
0.0027,
0.0025,
0.0026,
0.0032,
0.0063
]

############################################

W2013A = [180,
182,
183,
184,
186
]

W2013R = [5.3491,
5.3559,
5.3611,
5.3658,
5.3743
]

W2013dR = [0.0022,
0.0017,
0.002,
0.0023,
0.0026
]

Re2013A = [185,
187
]

Re2013R = [5.3596,
5.3698
]

Re2013dR = [0.0172,
0.0173
]

Os2013A = [184,
186,
187,
188,
189,
190,
192
]

Os2013R = [5.3823,
5.3909,
5.3933,
5.3993,
5.4016,
5.4062,
5.4126
]

Os2013dR = [0.0022,
0.0017,
0.0018,
0.0011,
0.0012,
0.0013,
0.0015
]

Ir2013A = [182,
183,
184,
185,
186,
187,
188,
189,
191,
193
]

Ir2013R = [5.3705,
5.378,
5.3805,
5.3854,
5.39,
5.3812,
5.3838,
5.3898,
5.3968,
5.4032
]

Ir2013dR = [0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061,
0.1061
]


Pt2013A = [183,
184,
185,
186,
187,
188,
189,
190,
191,
192,
193,
194,
195,
196,
198
]

Pt2013R = [5.4038,
5.4015,
5.4148,
5.4037,
5.4063,
5.4053,
5.406,
5.4108,
5.4102,
5.4169,
5.4191,
5.4236,
5.427,
5.4307,
5.4383
]

Pt2013dR = [0.0036,
0.0036,
0.0028,
0.0036,
0.0037,
0.0034,
0.0035,
0.003,
0.0031,
0.0028,
0.0027,
0.0025,
0.0026,
0.0027,
0.0032
]

IrEBIT2020A = [191,193
]

# IrEBIT2020R = [5.4422,
# 5.4486
# ]

IrEBIT2020dR = [0.0077,
    0.0076]

IrEBIT2020R = [5.4307,
    5.4371]

# IrEBIT2020R = [5.4307]
# IrEBIT2020dR = [0.0077]

plt.figure(figsize=(8.5,6)) 
plt.ylabel('RMS Charge radius [fm]')
plt.xlabel('Nucleon Number')
plt.minorticks_on()

#W2004A = [x-0.1 for x in W2004A]
plt.scatter(W2004A, W2004R, c='g', alpha=0.1, label='W (2004)')
plt.errorbar(W2004A, W2004R, yerr=W2004dR, ls='none', c='g', alpha=0.1, capsize=4)
plt.scatter(W2013A, W2013R, c='g', alpha=0.4, label='W (2013)')
plt.errorbar(W2013A, W2013R, yerr=W2013dR, ls='none', c='g', alpha=0.4 ,capsize=4)

#Re2004A2 = [x-0.15 for x in Re2004A]
plt.scatter(Re2004A, Re2004R, c='tab:orange', alpha=0.3, label='Re (2004)')
plt.errorbar(Re2004A, Re2004R, yerr=Re2004dR, ls='none', c='tab:orange', alpha=0.3,capsize=4)
plt.scatter(Re2013A, Re2013R, c='tab:orange', alpha=0.5, label='Re (2013)')
plt.errorbar(Re2013A, Re2013R, yerr=Re2013dR, ls='none', c='tab:orange', alpha=0.5,capsize=4)

#Os2004A = [x-0.2 for x in Os2004A]
plt.scatter(Os2004A, Os2004R, c='b', alpha=0.4, label='Os (2004)')
plt.errorbar(Os2004A, Os2004R, yerr=Os2004dR, ls='none', c='b', alpha=0.4,capsize=4)
plt.scatter(Os2013A, Os2013R, c='b', alpha=0.75, label='Os (2013)')
plt.errorbar(Os2013A, Os2013R, yerr=Os2013dR, ls='none', c='b', alpha=0.75,capsize=4)

#Ir2004A = [x-0.25 for x in Ir2004A]
plt.scatter(Ir2004A, Ir2004R,c='r', alpha=0.3, label='Ir (2004)')
plt.errorbar(Ir2004A, Ir2004R, yerr=Ir2004dR, c='r',ls='none', alpha=0.3,capsize=4)
# plt.scatter(Ir2013A, Ir2013R,c='r', alpha=0.6, label='Ir (2013)')
# plt.errorbar(Ir2013A, Ir2013R, yerr=Ir2013dR,c='r', ls='none', alpha=0.6,capsize=4)
plt.scatter(IrEBIT2020A, IrEBIT2020R, c='r', alpha=1, label='Ir (This work)')
plt.errorbar(IrEBIT2020A, IrEBIT2020R, yerr=IrEBIT2020dR, c='r', ls='none', alpha=1,capsize=4)

#Pt2004A = [x-0.3 for x in Pt2004A]
plt.scatter(Pt2004A, Pt2004R, c='tab:purple', alpha=0.3, label='Pt (2004)')
plt.errorbar(Pt2004A, Pt2004R, yerr=Pt2004dR, ls='none', c='tab:purple', alpha=0.3,capsize=4)
plt.scatter(Pt2013A, Pt2013R, c='tab:purple', label='Pt (2013)', alpha=0.7)
plt.errorbar(Pt2013A, Pt2013R, yerr=Pt2013dR, ls='none', c='tab:purple', alpha=0.7,capsize=4)

# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.legend( loc='lower right')
#plt.legend()
plt.xticks(np.arange(start=180, stop=200, step=2))
#plt.show()
plt.close()

#### deWitt MXS measurements and shifted values 


W2004DWA = [183]
#W2004DWA = [182.8]
W2004DWR = [5.3300]
W2004DWdR = [0.1500]

Os2004DWA = [189]
#Os2004DWA = [188.8]
Os2004DWR = [5.3600]
Os2004DWdR = [0.1500]

Ir2004DWA = [191,193]
#Ir2004DWA = [190.8,192.8]
Ir2004DWR = [5.3981, 5.4019]
Ir2004DWdR = [0.1061,0.1061]

W2013DWA = [183]
W2013DWR = [5.3611]
W2013DWdR = [0.0020]

Os2013DWA = [189]
Os2013DWR = [5.4016]
Os2013DWdR = [0.0012]

Ir2013DWA = [191, 193]
Ir2013DWR = [5.3968, 5.4032]
Ir2013DWdR = [0.1061, 0.1061]

# IrEBIT2020An = [191.2, 193.2]
IrEBIT2020An = [191, 193]
#IrEBIT2020An = [191]
csize = 4
# plt.figure() 
# plt.ylabel('RMS Charge radius [fm]')
# plt.xlabel('Nucleon Number')
# plt.minorticks_on()
# plt.scatter(W2004DWA, W2004DWR, label='$^{183}W$ (Angeli 2004)')
# plt.errorbar(W2004DWA, W2004DWR, yerr=W2004DWdR,capsize=csize, ls='none')

# plt.scatter(W2013DWA, W2013DWR, label='$^{183}W$ (Angeli 2013)')
# plt.errorbar(W2013DWA, W2013DWR, yerr=W2013DWdR,capsize=csize,  ls='none')

# plt.scatter(Os2004DWA, Os2004DWR, label='$^{189}Os$ (Angeli 2004)')
# plt.errorbar(Os2004DWA, Os2004DWR, yerr=Os2004DWdR,capsize=csize,  ls='none')


# plt.scatter(Os2013DWA, Os2013DWR, label='$^{189}Os$ (Angeli 2013)')
# plt.errorbar(Os2013DWA, Os2013DWR, yerr=Os2013DWdR,capsize=csize, ls='none')

# plt.scatter(Ir2004DWA, Ir2004DWR, label='$^{191,193}Ir$ (Angeli 2004)')
# plt.errorbar(Ir2004DWA, Ir2004DWR, yerr=Ir2004DWdR, capsize=csize, ls='none')


# plt.scatter(Ir2013DWA, Ir2013DWR, label='$^{191,193}Ir$ (Angeli 2013)')
# plt.errorbar(Ir2013DWA, Ir2013DWR, yerr=Ir2013DWdR, capsize=csize, ls='none')


# plt.scatter(Ir2013DWA, IrEBIT2020R, c='r', label='$^{191,193}Ir$ (This work)')
# plt.errorbar(Ir2013DWA, IrEBIT2020R, yerr=IrEBIT2020dR, c='r', capsize=csize,ls='none')


# plt.legend(loc='upper left', bbox_to_anchor=(0.1, 0.8))
# #plt.show()
# plt.close


fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
fig.subplots_adjust(hspace=0.005, wspace= 0.03)

ax1.scatter(W2004DWA, W2004DWR, c='g', marker="x", alpha = 0.25, label='$^{183}W$ (Angeli 2004)')
ax1.errorbar(W2004DWA, W2004DWR, c='g', alpha = 0.25,yerr=W2004DWdR,capsize=csize, ls='none')
ax1.scatter(W2013DWA, W2013DWR, c='g', label='$^{183}W$ (Angeli 2013)')
ax1.errorbar(W2013DWA, W2013DWR, c='g', yerr=W2013DWdR,capsize=csize,  ls='none')

ax2.scatter(W2004DWA, W2004DWR, c='g', marker="x",alpha = 0.25,label='$^{183}W$ (Angeli 2004)')
ax2.errorbar(W2004DWA, W2004DWR, c='g', alpha = 0.25,yerr=W2004DWdR,capsize=csize, ls='none')
ax2.scatter(W2013DWA, W2013DWR, c='g',label='$^{183}W$ (Angeli 2013)')
ax2.errorbar(W2013DWA, W2013DWR, c='g',yerr=W2013DWdR,capsize=csize,  ls='none')


ax2.scatter(Os2004DWA, Os2004DWR, c='b', marker="x",alpha=0.25, label='$^{189}Os$ (Angeli 2004)')
ax2.errorbar(Os2004DWA, Os2004DWR, c='b', alpha=0.25,yerr=Os2004DWdR,capsize=csize,  ls='none')

ax2.scatter(Os2013DWA, Os2013DWR, c='b', label='$^{189}Os$ (Angeli 2013)')
ax2.errorbar(Os2013DWA, Os2013DWR, c='b', yerr=Os2013DWdR,capsize=csize, ls='none')

ax2.scatter(Ir2004DWA, Ir2004DWR, c='r', marker="x",alpha=0.25, label='$^{191,193}Ir$ (Angeli 2004)')
ax2.errorbar(Ir2004DWA, Ir2004DWR, c='r', alpha = 0.25, yerr=Ir2004DWdR, capsize=csize, ls='none')

# ax2.scatter(Ir2013DWA, Ir2013DWR, label='$^{191,193}Ir$ (Angeli 2013)')
# ax2.errorbar(Ir2013DWA, Ir2013DWR, yerr=Ir2013DWdR, capsize=csize, ls='none')

ax2.scatter(IrEBIT2020An, IrEBIT2020R, c='r', label='$^{191,193}Ir$ (This work)')
ax2.errorbar(IrEBIT2020An, IrEBIT2020R, yerr=IrEBIT2020dR, c='r', capsize=csize,ls='none')

ax1.set_xlim(182,184)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xticklabels(labels=[str('\n') ,183, str('\n')], rotation=0)
ax2.set_xlim(188,194)
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.set_xticklabels(labels=[str('\n') ,189, str('\n'), 191, str('\n'), 193, str('\n')], rotation=0)
#ax2.set_xticklabels(labels=[str('\n') ,189, str('\n'), 191, str('\n')], rotation=0)

ax1.spines.right.set_visible(False)
ax2.spines.left.set_visible(False)
ax1.yaxis.tick_left()
ax1.tick_params(labelright='off')
ax1.yaxis.tick_left()
ax2.tick_params(labelleft='off')
ax2.yaxis.tick_right()
ax1.minorticks_on()
ax2.minorticks_on()
#ax2.yaxis('off')
# d = .5  # proportion of vertical to horizontal extent of the slanted line
# kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
#               linestyle="none", color='k', mec='k', mew=1, clip_on=False)
# ax1.plot([1, 0], transform=ax1.transAxes, **kwargs)
# ax2.plot([0, 0], transform=ax2.transAxes, **kwargs)
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
ax2.plot((-d, +d), (-d, +d), **kwargs)
ax1.set_ylabel('Nuclear radius [fm]')
fig.supxlabel('Nucleon number')
#plt.show()
plt.close()


methNum = [1,2,3,4,5,6,7,8]
methlab = ['A', 'B', 'C', 'D', 'F', 'G', 'H', 'I']
measSep = [3.501, 3.502, 3.508, 3.506, 3.497, 3.502, 3.495, 3.499]
dmeasSep = [0.007, 0.009, 0.007, 0.013, 0.005, 0.022, 0.007, 0.021]

R191 = [5.4339,
    5.436,
    5.423,
    5.4274,
    5.4252,
    5.4143,
    5.4295,
    5.4208]
dR191 = [0.006722351,
    0.009991997,
    0.010220342,
    0.017960484,
    0.007239475,
    0.02903739,
    0.012492758,
    0.032636942]
plt.figure() 
plt.scatter(methlab, R191, s=6)
plt.errorbar(methlab, R191, dR191, capsize=4, ls='none')

plt.scatter(3.5, 5.4339, c='g')
plt.errorbar(3.5, 5.4339, 0.0101, c='g', capsize=4, ls='none')

plt.scatter(4.5, 5.4307, c='r')
plt.errorbar(4.5, 5.4307, 0.0077, c='r', capsize=4, ls='none')


plt.minorticks_on()
plt.xlabel('Method')
plt.ylabel('Nuclear charge radius [fm]')
plt.show()
plt.close()

dat = pd.read_csv(r'C:\\Users\\ahosi\\Downloads\\SpectraData.csv')
plt.figure(figsize=(20,7))
plt.plot(dat['Wavelength (nm)'], dat['Os_sum'])
plt.plot(dat['Wavelength (nm)'],dat['Ir_sum'])
plt.minorticks_on()
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (arb)')
plt.xlim(left=4, right=19.5)
plt.ylim(top=19000, bottom=8330)
#plt.show()
plt.close()


fig, ax1 = plt.subplots(figsize=(15,7))
ax1.plot(dat['Wavelength (nm)'], dat['Os_sum'])
ax1.plot(dat['Wavelength (nm)'], dat['Ir_sum'])
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (arb)')
ax1.minorticks_on()
ax1.set_xlim(left=4, right=19.5)
ax1.set_ylim(bottom=8000)
ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.35, 0.35, 0.3, 0.5])
ax2.set_axes_locator(ip)

ax2.plot(dat['Wavelength (nm)'], dat['Os_sum'])
ax2.plot(dat['Wavelength (nm)'], dat['Ir_sum'])
mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')
ax2.set_xlim(left=7.2, right=7.8)
ax2.minorticks_on()
ax2.set_ylim(bottom=8500, top=14000)

plt.show()