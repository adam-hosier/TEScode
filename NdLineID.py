import mass 
import numpy as np 
import mass.materials
from uncertainties import unumpy as unp
from uncertainties import ufloat
import pylab as plt 
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os 
#import ebit_util
import pandas as pd
import matplotlib.colors as mcolors
from scipy.stats import binned_statistic_2d
from scipy import stats
from fit_utils import MultiPeakGaussian
from itertools import cycle 
import glob 
import os 
from lmfit import minimize, Parameters, report_fit, Model, Parameter
from lmfit.models import GaussianModel
from lmfit.models import SplitLorentzianModel 
from lmfit.models import VoigtModel 
from lmfit.models import PseudoVoigtModel
from lmfit.models import SkewedVoigtModel 


#### Using actual theory convoluted spectra generated by Yuri 
tEl = 'Nd_'
cxfactor = '1e12'
ftheoryloc = str('C:\\data\\theory')
ftname = 'conv1.4_SP1000.00_1.00e+12_0_0_2000.dat'
tdat3 = pd.read_fwf(ftheoryloc+'\\'+str(tEl)+cxfactor+'\\'+str(ftname))
tdat3 = np.array(tdat3)

ftname4 = 'conv1.4_SP1000.00_1.00e+10_0_0_2000.dat'
tdat4 = pd.read_fwf(ftheoryloc+'\\'+str(tEl)+cxfactor+'\\'+str(ftname4))
tdat4 = np.array(tdat4)

ftname5 = 'conv1.4_SP1000.00_1.00e+11_0_0_2000.dat'
tdat5 = pd.read_fwf(ftheoryloc+'\\'+str(tEl)+cxfactor+'\\'+str(ftname5))
tdat5 = np.array(tdat5)

ftname6 = 'conv1.4_SP1000.00_2.00e+10_0_0_2000.dat'
tdat6 = pd.read_fwf(ftheoryloc+'\\'+str(tEl)+cxfactor+'\\'+str(ftname6))
tdat6 = np.array(tdat6)


iffit = False



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#floc = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022')
floc = str('C:\\data\\TES_Spectra_1eVbin')
floc = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\TES_Calibration_Lines\\20221221')
floc2 = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\TES_Calibration_Lines\\SDr_Test\\')
floc3 = str('C:\\data\\TES_newcal')

f3names = glob.glob(floc3+'\\*.csv')


df3 = dict() 
for fpath in f3names: 
    df3[os.path.basename(fpath)[:-4]] = pd.read_csv(r""+fpath)


ftheoryloc = str('C:\\data\\theory')
teseffloc = str('C:\\data')
theoryEl = 'Nd'
eltitle = theoryEl
#floc = str('C:\\data\\calibratedTES_Dec2022')
ddest = str('C:\\data\\Line_ID_Nd')
date = str('202212')
day = str('19')
runnum = str('0000')

#statelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']
#statelist = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]
statelist = ['H', 'K', 'P', 'W', 'Y', 'AA']

coAdd = False
df = dict()

Cnames = ['energy', 'Intensity', 'Spec Charge', 'con1i', 'con2i', 'Numberi', 'Ji', '-',
          'con1f', 'con2f', 'Numberf', 'Jf', ':', 'Intensity2', '|']
tdat = pd.read_csv(r""+ftheoryloc+'\\'+theoryEl+'.csv', names=Cnames)

#newcal = pd.read_csv(r"C:\\data\\TES_ReCaltest_Calibration.csv")


tenergy = tdat['energy']
tintensity = tdat['Intensity']

  
arrx = df3['20221219_0000_AC_cal']['20221219_0000_AA_Energy']
arry = df3['20221219_0000_AC_cal']['20221219_0000_AA_Counts']


########## Fitting of lines for Line ID individually 

def datcoll(y, x, c, r):  #in single spectra loop
    res = dict()

    ydat = []   
    ydate = []                   
    wavedat = []
    pix = []
    bg = 0                  #background outside of radius of data 
    photc = []                  #photon count

    #Collecting data around a point of interest defined by 'c' and 'r' above, along with centroid calcs
    for i in range(1,2*r):
        ydat.append(y[c - r + i])
        wavedat.append(x[c-r+i])
        pix.append(c-r+i)
   
    ydat = np.array(ydat)
    ydate = np.sqrt(ydat)
    wavedat = np.array(wavedat)
    pix = np.array(pix) 

    res['y'] = ydat
    res['ye'] = ydate 
    res['x'] = pix 
    res['wave'] = wavedat

    return res 


def vfit(x, y, E, r, num_peaks=1): 
    rez = dict() 
    ydat = []
    xdat = []
    c = find_nearest(x, E)
    for k in range(1,2*r):
        ydat.append(y[c - r + k])
        xdat.append(x[c-r+k])

    ydat = np.array(ydat)
    xdat = np.array(xdat)
    pars = Parameters()
    for i in range(num_peaks):
        if i == 0 and num_peaks > 1:
            rez['Composite model'] = VoigtModel(prefix='V'+str(i+1)+'_')
        elif i>0 and num_peaks > 1: 
            rez['Composite model'] += VoigtModel(prefix='V'+str(i+1)+'_')
        
        rez['Voigt_mod_'+str(i+1)] = VoigtModel(prefix='V'+str(i+1)+'_') 

        pars1 = rez['Voigt_mod_'+str(i+1)].make_params()
        pars1['V'+str(i+1)+'_center'].set(min=0.85*np.min(xdat), max=1.15*np.max(xdat), value = np.mean(xdat), vary=True)
        pars1['V'+str(i+1)+'_fwhm'].set(min=4.715, max=5.50, value=4.9475, vary=False)
        #pars1['V'+str(i+1)+'_fwhm'].set(min=4.715, max=5.50, vary=True)
        # if i>0 and num_peaks>1: 
        #     pars1['V'+str(i+1)+'_sigma'].set(expr='V1_sigma')

        pars.update(pars1)

    if num_peaks > 1: 
        modtemp = rez['Composite model']
    elif num_peaks == 1:
        modtemp = rez['Voigt_mod_1']
    
    out = modtemp.fit(ydat, pars, x=xdat, weights=np.sqrt(ydat), nan_policy='omit')
    pars.update(out.params)

    neweval = modtemp.eval(pars, x=np.linspace(np.min(xdat), np.max(xdat), num=1000))
    comps = out.eval_components(x=np.linspace(np.min(xdat), np.max(xdat), num=1000))
    rez['Params'] = pars 
    rez['out'] = out
    rez['neweval'] = neweval 
    rez['xdat'] = xdat 
    rez['ydat'] = ydat 
    rez['num_data_points'] = k 
    rez['newevalx'] = np.linspace(np.min(xdat), np.max(xdat), num=1000)

    return rez 

plottitle = '20221221_0002_AIOAHcal_AF_Counts'
xd = df3['20221221_0002_AIOAH_cal']['20221221_0002_A_Energy']
yd = df3['20221221_0002_AIOAH_cal'][plottitle]


# test1 = vfit(xd, yd, 1204, 20, num_peaks=1)
test2 = vfit(xd, yd, 1203, 15, num_peaks=1)
# test3 = vfit(xd, yd, 1172, 10, num_peaks=1)
#test4 = vfit(xd, yd, 1361, 7, num_peaks=1)
# test5 = vfit(xd, yd, 843, 12, num_peaks=2)
# test6 = vfit(xd, yd, 1522, 10, num_peaks=1)
# test7 = vfit(xd, yd, 1546, 10, num_peaks=1)
#print(test1['Params'].pretty_print())
#print(test1['out'].fit_report())

# print('######################')
# print(test2['Params'].pretty_print())
# print('######################')
# print(test3['Params'].pretty_print())
# print('######################')
# print(test4['Params'].pretty_print())
# print('######################')
# print(test5['Params'].pretty_print())
# print('######################')
# print(test6['Params'].pretty_print())
# print('######################')
# print(test7['Params'].pretty_print())
# print('######################')


print(test2['out'].fit_report())

plt.figure() 
plt.plot(xd, yd, c='r', label='data')
# plt.plot(test1['newevalx'], test1['neweval'], c='b', label='fit')
# plt.plot(test2['newevalx'], test2['neweval'], c='b')
# plt.plot(test3['newevalx'], test3['neweval'], c='b')
plt.plot(test2['newevalx'], test2['neweval'], c='b')
# plt.plot(test5['newevalx'], test5['neweval'], c='b')
# plt.plot(test6['newevalx'], test6['neweval'], c='b')
# plt.plot(test7['newevalx'], test7['neweval'], c='b')
plt.xlabel('Energy (eV)')
plt.ylabel('Counts per 1 eV bin')
plt.title(str(plottitle))
plt.legend()
plt.show()
plt.close() 
