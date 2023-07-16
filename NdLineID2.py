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
from lmfit.models import LinearModel 


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


def vfit(x, y, E, r, num_peaks=1, linbackground = False): 
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
        pars1['V'+str(i+1)+'_center'].set(min=np.min(xdat), max=np.max(xdat), value = E, vary=True)
        #pars1['V'+str(i+1)+'_fwhm'].set(min=4.715, max=5.50, value=4.9475, vary=False)
        pars1['V'+str(i+1)+'_amplitude'].set(vary=True, min=0)
        #pars1['V'+str(i+1)+'_height'].set(vary=True, min=0)
        pars1['V'+str(i+1)+'_fwhm'].set(min=4)
        pars1['V'+str(i+1)+'_sigma'].set(value=2, min=0,  vary=True)
        #pars1['V'+str(i+1)+'_sigma'].set(min=1.3, max=1.50, value=1.9, vary=False)
        pars1['V'+str(i+1)+'_gamma'].set(value = 0, min =0, vary=False)

        if i>0 and num_peaks>1: 
            pars1['V'+str(i+1)+'_sigma'].set(expr='V1_sigma')
            #pars1['V'+str(i+1)+'_gamma'].set(expr='V1_gamma')

            pars1['V'+str(i+1)+'_center'].set(expr ='V1_center+6.1', vary=True)
        

        pars.update(pars1)


    if num_peaks > 1: 
        if linbackground is True: 
            rez['Composite model'] += LinearModel(prefix='Background_')
            pars2 = LinearModel(prefix='Background_').make_params()
            pars2['Background_slope'].set(value=0, vary=True)
            pars2['Background_intercept'].set(value=np.min(ydat), vary=True)
            pars.update(pars2)
        modtemp = rez['Composite model']
    elif num_peaks == 1:
        if linbackground is True: 
            rez['Voigt_mod_1'] += LinearModel(prefix='Background_')
            pars2 = LinearModel(prefix='Background_').make_params()
            pars2['Background_slope'].set(value=0, vary=True)
            pars2['Background_intercept'].set(value=0, vary=True)
            pars.update(pars2)
        modtemp = rez['Voigt_mod_1']
    
    #out = modtemp.fit(ydat, pars, x=xdat, weights=np.sqrt(ydat), nan_policy='omit')
    out = modtemp.fit(ydat, pars, x=xdat, weights=np.sqrt(ydat))
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
    rez['comps'] = comps

    return rez 

##energy dependent 
# expstatelist = ['H', 'K', 'P', 'W', 'Y', 'AA']
# plottitle = '20221219_0000_Y_Counts'
# xd = df3['20221219_0000_BCDAC_cal']['20221219_0000_B_Energy']
# yd = df3['20221219_0000_BCDAC_cal'][plottitle]



#### T V X Z AB AD AF
##Density dependence data
#expstatelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']
expstatelist = ['T', 'AF', 'AD', 'AB', 'V', 'Z', 'X']
expstatelist2 = ['H', 'K', 'P', 'W', 'Y', 'AA']
plottitle = '20221221_0002_AIOAHcal_T_Counts'
xd = df3['20221221_0002_AIOAH_cal']['20221221_0002_A_Energy']
yd = df3['20221221_0002_AIOAH_cal'][plottitle]
# beamcurr = [36.80,
#     18.40,
#     9.20,
#     13.80,
#     23.00,
#     27.60,
#     32.20]
beamcurr = [36.80,
    32.20,
    27.60,
    23.00,
    18.40,
    13.80,
    9.20]

t_norm = [3163, 
    1769, 
    1394, 
    650, 
    328, 
    365, 
    3109]

npeak = 2
lenergy = 1201
npoints = 5
lback = False
test2 = vfit(xd, yd, lenergy, npoints, num_peaks=npeak, linbackground=lback)
rescolnames = ['line number', 'center', 'center unc', 'height', 'height unc', 'FWHM', 'FWHM unc', 'sigma', 'sigma unc', 'gamma', 'gamma unc']
touttest = np.zeros((1,11))
for state in expstatelist:
    ptitle = str('20221221_0002_AIOAHcal_')+str(state)+str('_Counts')
    xd = df3['20221221_0002_AIOAH_cal']['20221221_0002_A_Energy']
    yd = df3['20221221_0002_AIOAH_cal'][ptitle]
    test3 = vfit(xd, yd, lenergy, npoints, num_peaks=npeak, linbackground=lback)
    paramlist = test3['Params'] 
    outtest = np.zeros((1,11))
    for ln in range(1,npeak+1): 
        outdat = np.array([[ln, 
                          paramlist['V'+str(ln)+'_center'].value, paramlist['V'+str(ln)+'_center'].stderr, 
                          paramlist['V'+str(ln)+'_height'].value, paramlist['V'+str(ln)+'_height'].stderr, 
                          paramlist['V'+str(ln)+'_fwhm'].value, paramlist['V'+str(ln)+'_fwhm'].stderr,
                          paramlist['V'+str(ln)+'_sigma'].value, paramlist['V'+str(ln)+'_sigma'].stderr, 
                          paramlist['V'+str(ln)+'_gamma'].value, paramlist['V'+str(ln)+'_gamma'].stderr]])
        
        touttest = np.vstack((touttest, outdat))
    

prepdata = pd.DataFrame(data = touttest, columns=rescolnames)
prepdata.to_excel('C:\\Users\\ahosi\\anaconda3\\envs\\TESenv\\TEScode\\testout.xlsx', index=False)

plottitle = '20221221_0002_AIOAHcal_T_Counts'
xd = df3['20221221_0002_AIOAH_cal']['20221221_0002_A_Energy']
yd = df3['20221221_0002_AIOAH_cal'][plottitle]


for state in expstatelist2: 
    ptitle = str('20221219_0000_')+str(state)+str('_Counts')




plt.figure() 
tc = 0
for state in expstatelist: 
    ptitle = str('20221221_0002_AIOAHcal_')+str(state)+str('_Counts')
    xd = df3['20221221_0002_AIOAH_cal']['20221221_0002_A_Energy']
    #yd = df3['20221221_0002_AIOAH_cal'][ptitle] / t_norm[tc]
    yd = df3['20221221_0002_AIOAH_cal'][ptitle] / np.max(df3['20221221_0002_AIOAH_cal'][ptitle])
    #yd = df3['20221221_0002_AIOAH_cal'][ptitle] 
    plt.plot(xd, yd, label='2.04 keV, '+str(beamcurr[tc])+' mA')
    tc +=1 
#plt.plot(xd, yd, c='r', label='data')
# plt.plot(test1['newevalx'], test1['neweval'], c='b', label='fit')
# plt.plot(test2['newevalx'], test2['neweval'], c='b')
# plt.plot(test3['newevalx'], test3['neweval'], c='b')

#plt.plot(test2['newevalx'], test2['neweval'], c='b')

# plt.plot(test5['newevalx'], test5['neweval'], c='b')
# plt.plot(test6['newevalx'], test6['neweval'], c='b')
# plt.plot(test7['newevalx'], test7['neweval'], c='b')

###plotting theory 
#plt.plot(tdat3[:,0], np.max(yd)*tdat3[:,1]/np.max(tdat3[:,1]), label=ftname)
# plt.plot(tdat4[:,0], np.max(yd)*tdat4[:,1]/np.max(tdat4[:,1]), label=ftname4)
# plt.plot(tdat5[:,0], np.max(yd)*tdat5[:,1]/np.max(tdat5[:,1]), label=ftname5)
# plt.plot(tdat6[:,0], np.max(yd)*tdat6[:,1]/np.max(tdat6[:,1]), label=ftname6)
# qlines = [1204, 1545.7, 1712, 1828, 1899, 1948, 1983]
# qlines2 = [1240, 1606, 1733, 1876, 1953, 2007, 2047]
# for t in qlines: 
#     plt.axvline(x=t, c='k', ls='--')

# for t in qlines2: 
#     plt.axvline(x=t, c='g', ls='--')


fitted_lines = [842.64,
    865.94,
    876.00,
    884.12,
    932.62,
    942.75,
    964.29,
    976.96,
    984.83,
    999.39,
    1151.03,
    1157.89,
    1172.29,
    1203.89,
    1241.34,
    1360.30,
    1379.62,
    1386.26,
    1395.85,
    1461.77,
    1475.76,
    1494.18,
    1502.47,
    1522.32,
    1546.00,
    1606.39,
    1711.60,
    1733.94,
    1825.05,
    1898.39,
    1984.78,
    1948.14]

tbd_lines = [1415.83,
    1424.17,
    1557.25,
    1572.99,
    1580.66,
    1593.13,
    1652.47,
    1774.46,
    1781.54,
    1807.04,
    1847.57,
    1920.77,
    1929.08,
    1972.84]


lcount = 1

yd = np.array(yd)
# for l in fitted_lines: 
#     #plt.axvline(x=l, c='k', ls='--')
#     # plt.text(l, yd[find_nearest(xd, l)]+300, str(lcount))
#     #plt.text(l, yd[find_nearest(xd, l)]+0.02, str(lcount))
#     plt.text(l, 1.025+0.025*(-1)**(lcount), str(lcount))
#     lcount += 1 

# for l in tbd_lines:
#     # plt.text(l, yd[find_nearest(xd, l)]+300, str(lcount), c='g')
#     plt.text(l, yd[find_nearest(xd, l)], str(lcount), c='g')
#     lcount += 1

# for l in range(npeak):
#     plt.plot(test2['newevalx'], test2['comps']['V'+str(l+1)+'_'], label='Voigt #'+str(l+1))

plt.xlabel('Energy (eV)')
plt.ylabel('Normalized Intensity')
#plt.title(str(plottitle))
#plt.axhline(y=0, c='k', ls='--')
plt.ylim(bottom=0)
plt.legend()
plt.minorticks_on()
# plt.xlim((np.min(test2['newevalx']), np.max(test2['newevalx'])))
# plt.ylim((0, 1.05*np.max(test2['ydat'])))
#plt.xlim((1125, 1300))
plt.xlim((700, 2100))
plt.show()
plt.close() 
#print(plottitle)