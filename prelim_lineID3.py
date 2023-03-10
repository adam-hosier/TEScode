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

theoryres = np.shape(tdat3)[0]
minE = np.min(tdat3[:,0])
maxE = np.max(tdat3[:,0])
resStep = (maxE-minE)/theoryres
xray_range = np.arange(minE, maxE, resStep)
EBIT_model_old = mass.materials.filterstack_models['EBIT 2018']
EBIT_model_new = mass.materials.FilterStack(name='2022 Stack')
EBIT_model_new.add_Film(name='Electroplated Au Absorber', material='Au',area_density_g_per_cm2=ufloat(0.00186,0.00006), fill_fraction=ufloat(1.000,0), absorber=True)
EBIT_model_new.add_AlFilmWithOxide(name='50mK Filter', Al_thickness_nm=100)
EBIT_model_new.add_AlFilmWithOxide(name='5K Filter', Al_thickness_nm=112)
EBIT_model_new.add_AlFilmWithOxide(name='50K Filter', Al_thickness_nm=116)
#EBIT_model_new.add_Film(name='Ni Mesh', material='Ni', area_density_g_per_cm2=ufloat(0.0134,0.0018), fill_fraction=ufloat(0.170,0.010), absorber=False)
EBIT_model_new.add_LEX_HT(name='Luxel Window #1')
EBIT_model_new.add_LEX_HT(name='Luxel Window #2')

eff1 = EBIT_model_new.get_efficiency(xray_range)
eff_unc1 = EBIT_model_new.get_efficiency(xray_range, uncertain=True)
uncs1 = unp.std_devs(eff_unc1)
#EBIT_model_new.plot_efficiency(xray_energies_eV=xray_range)

effdat1 = np.vstack((xray_range, eff1, uncs1)).T

dfeff = pd.DataFrame(data=effdat1, columns=['Energy (eV)', 'Efficiency %', 'Efficiency % Uncertainty'])


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
beamen = [2.04,
    2.04,
    2.04,
    2.04,
    2.04,
    2.04,
    2.04]

beamcurr = [36.8,
    18.4,
    9.2,
    13.8,
    23,
    27.6,
    32.2]
coAdd = False
df = dict()
# minenergy = 500
# maxenergy = 5000
# binsize = 1
# numbins = int(np.round((maxenergy-minenergy)/binsize))
Cnames = ['energy', 'Intensity', 'Spec Charge', 'con1i', 'con2i', 'Numberi', 'Ji', '-',
          'con1f', 'con2f', 'Numberf', 'Jf', ':', 'Intensity2', '|']
tdat = pd.read_csv(r""+ftheoryloc+'\\'+theoryEl+'.csv', names=Cnames)

#newcal = pd.read_csv(r"C:\\data\\TES_ReCaltest_Calibration.csv")

#sDrdat = pd.read_csv(r"C:\\data\\theory\\SDR.csv")
#sDrdat = pd.read_csv(r"C:\\Users\\ahosi\Downloads\\conv_1.00000000E-15.csv")
#newdat2 = pd.read_csv(r""+floc2+'AlSiClK_only.csv')
#newdat3 = pd.read_csv(r""+floc2+'AlSiTiFe_only.csv')
#effdat = pd.read_csv(r""+teseffloc+'\\'+'TES_Efficiency_Dec2022.csv')
#teseff = effdat['Efficiency %']
teseff = dfeff['Efficiency %']
tenergy = tdat['energy']
tintensity = tdat['Intensity']


# for s in statelist: 
#     state = str(s)
#     df[state] = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'_sum'+'.csv')
#     #df = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'.csv')

#     # counts, bin_edges = np.histogram(df[state]['energy'], bins=numbins, range=(minenergy, maxenergy))
#     # df[state+str(' counts')]= counts
#     # df[state+str(' bin_edges')] = bin_edges


# energy = df['energy']
# time = df['time']

# counts, bin_edges = np.histogram(energy, bins=numbins, range=(minenergy, maxenergy))

###################################
# plt.figure()
# plt.ylabel('Counts per '+str(binsize)+' eV bin')
# #plt.ylim(bottom=0, top=1.1*np.max(counts))
# #plt.xlim(left=3075, right=3175)
# plt.xlabel('energy (eV)')
# plt.title(date+day+'_'+runnum)
# #plt.title(date+day+'_'+runnum+'_'+state)
# for state in statelist: 
#     plt.plot(df[state+str(' bin_edges')][:-1], df[state+str(' counts')], label=state)
# plt.legend()
# #plt.show() 
# ##################################
#tsig = 1.0333       #standard deviation of theoretical gaussian assuming 4.5 eV FWHM instrument resolution

state = 'AA'
# arry = df[state+str(' counts')]
# arrx = df[state+str(' bin_edges')][:-1]

# arry = df[state]['1']
# arrx = df[state]['0']

# arry = df[state]['Counts']
# arrx = df[state]['Energy [eV]']   
arrx = df3['20221219_0000_AC_cal']['20221219_0000_AA_Energy']
arry = df3['20221219_0000_AC_cal']['20221219_0000_AA_Counts']

tsig = 1.877302      #stand dev of theo gauss assuming ~6.7 eV FWHM instrument resolution 
def tgauss(x, A, mu, sig):
    return A*np.exp(-(x-mu)**2 / (2*sig**2))

tplot = np.linspace(500, 8000, num=theoryres)
#norm = np.max(tintensity)   
#tintensity = tintensity / norm


##Normalizing over all lines 
tx = tplot

ty = np.zeros((np.shape(tplot,)))
for i in range(np.shape(tenergy)[0]): 
    ty += tgauss(tplot, tintensity[i], tenergy[i], tsig)

ty = ty*teseff
norm = np.max(ty)
ty = ty / norm

### Normalizing per charge state 
cstates = tdat['Spec Charge'].unique() 
dfspec = dict() 
cnorm = []
normE = []
normcounts = []
dfspec['Tx'] = tx
m =0
for k in cstates:
    ctemp = np.zeros((0,np.shape(tdat)[1]))
    tdat2 = np.array(tdat)
    for i in range(np.shape(tdat)[0]):
        if tdat['Spec Charge'][i] == k: 
            ctemp = np.vstack((ctemp, tdat2[i,:]))
        else: 
            pass

    dfspec[str(k)] = ctemp 
    cnorm.append(np.max(ctemp[:,1]))
    normEindex = np.where(ctemp==np.max(ctemp[:,1]))[0][0] 
    normE.append(ctemp[normEindex, 0])
    if m == 2:
        expEindex = find_nearest(arrx, 1241)
    elif m==1: 
        expEindex = find_nearest(arrx, 1204)
    else:
        expEindex = find_nearest(arrx, ctemp[normEindex, 0])
    expEnorm = arry[expEindex]
    normcounts.append(expEnorm)
    ty2 = np.zeros((np.shape(tplot,)))
    for i in range(np.shape(ctemp)[0]): 
        ty2 += tgauss(tplot, ctemp[i,1], ctemp[i,0], tsig)

    ty2 = ty2*teseff
    norm2 = cnorm[m]  
    ty2 = normcounts[m]*ty2 / norm2 
    dfspec[str(k)+' calcint'] = ty2

    m += 1

lineSlist = ["-","--","-.",":"]

if iffit is True:
    res_dict = dict()
    a = MultiPeakGaussian(arr = arry, xs = arrx, num_peaks=25, num_poly=3)
    a.fit(return_dict = res_dict, same_sigma=True, function='voigt')
    t1 = res_dict['rez']
    #a.plot_fit(normalize_background=True)

    b = t1[1]
    c = np.zeros([1,9])
    for i in range(np.shape(t1)[0]):

        tempa = np.zeros([1,1])
        for l in range(np.shape(t1)[1]):
            temp = np.array([t1[i][l]])

            tempa = np.hstack((tempa,temp))
        
        c = np.vstack((c, tempa))

    c = c[1:,:]
    c = c[:,1:]

    lineIDdata = pd.DataFrame(data=c, columns=['center [eV]', 'center unc [eV]','height [counts]','height unc [counts]', 
        'std dev [eV]','std dev unc [eV]', 'FWHM [eV]', 'FWHM unc [eV]'])

    lineIDdata.to_csv(ddest+'/'+str(date)+str(day)+'_'+str(runnum)+'_'+str(state)+'.csv', index=True)


ebeamen = beamen[statelist.index(str(state))]
ebeamcurr = beamcurr[statelist.index(str(state))]
if iffit is True:
    xfit = res_dict['x0']
    yfit = res_dict['y0']
lcycler = cycle(lineSlist)
#statelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']

#statelist = ['T', 'V', 'X']
#statelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']
#statelist = ["E", "G", "K", "M", "Q", "R", "T", "V", "X", "Z", "AB", "AD", "AF"]
statelist = ['AA']



#statelist = ['AA']
plt.figure() 
#plt.title(str(date)+'_'+str(day)+'_'+str(runnum)+'_'+str(state)+' /// '+str(eltitle)+' '+ str(ebeamen)+' keV , '+str(ebeamcurr)+' mA')
plt.ylabel('Photon counts')
plt.xlabel('Energy (eV)')
plt.minorticks_on() 
#plt.plot(arrx, arry, c='b', label='Experimental data')
#norm1 = np.max(df[state]['Counts'])

# for state1 in statelist:
#     # ebeamen = beamen[statelist.index(str(state1))]
#     # ebeamcurr = beamcurr[statelist.index(str(state1))]
#     #plt.plot(df[state1]['0'], norm1*df[state1]['1']/np.max(df[state1]['1']), label=str(eltitle)+' '+ str(ebeamen)+' keV , '+str(ebeamcurr)+' mA')
#     plt.plot(df[state1]['Energy [eV]'], norm1*df[state1]['Counts']/np.max(df[state1]['Counts']), label=state1)
    
#     #plt.plot(df[state1]['0'], norm1*df[state1]['1']/np.max(df[state1]['1']), label=state1+' MASS only calibration')
#     # plt.plot(newcal['calibration'], norm1*df[state1]['1']/np.max(df[state1]['1']), label=state1+' + traditional calibration')
#     plt.plot(sDrdat['energy'], 505.8*sDrdat['intensity']/np.max(sDrdat['intensity']), label='S DR Theory')


if iffit is True:
    
    plt.plot(xfit, yfit, c='g', label='Fitted data')
    plt.plot(arrx, arry, c='b', label='data')
# plt.plot(tdat3[:,0], np.max(arry)*(tdat3[:,1]*eff1)/np.max(tdat3[:,1]*eff1),ls='--', label=ftname)
# plt.plot(tdat4[:,0], np.max(arry)*(tdat4[:,1]*eff1)/np.max(tdat4[:,1]*eff1),ls='--', label=ftname4)
# plt.plot(tdat5[:,0], np.max(arry)*(tdat5[:,1]*eff1)/np.max(tdat5[:,1]*eff1),ls='--', label=ftname5)
# plt.plot(tdat6[:,0], np.max(arry)*(tdat6[:,1]*eff1)/np.max(tdat6[:,1]*eff1),ls='--', label=ftname6)

# plt.plot(tx, ty*np.max(arry), c='r', label='Theory (all charge states)')
# for k in cstates:
#     plt.plot(dfspec['Tx'], dfspec[str(k)+' calcint'],c=np.random.rand(3,),ls=next(lcycler), label='Spectroscopic charge: '+str(k))
# plt.xlim((350, 1000))
# plt.ylim((0,500))
plt.legend() 
#plt.show() 
plt.close() 

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
        pars1['V'+str(i+1)+'_center'].set(min=0.85*np.min(xdat), max=1.15*np.max(xdat), value=np.median(xdat))
        if i>0 and num_peaks>1: 
            pars1['V'+str(i+1)+'_sigma'].set(expr='V1_sigma')

        pars.update(pars1)

    if num_peaks > 1: 
        modtemp = rez['Composite model']
    elif num_peaks == 1:
        modtemp = rez['Voigt_mod_1']
    
    out = modtemp.fit(ydat, pars, x=xdat)
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

plottitle = '20221221_0002_AIOAHcal_T_Counts'
xd = df3['20221221_0002_AIOAH_cal']['20221221_0002_A_Energy']
yd = df3['20221221_0002_AIOAH_cal']['20221221_0002_AIOAHcal_T_Counts']


test1 = vfit(xd, yd, 1204, 20, num_peaks=4)
test2 = vfit(xd, yd, 1241, 13, num_peaks=2)
test3 = vfit(xd, yd, 1172, 10, num_peaks=1)
test4 = vfit(xd, yd, 1360, 10, num_peaks=1)
test5 = vfit(xd, yd, 843, 12, num_peaks=2)
test6 = vfit(xd, yd, 1522, 10, num_peaks=1)
test7 = vfit(xd, yd, 1546, 10, num_peaks=1)
print(test1['Params'].pretty_print())
print(test1['out'].fit_report())

print('######################')
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

plt.figure() 
plt.plot(xd, yd, c='r', label='data')
plt.plot(test1['newevalx'], test1['neweval'], c='b', label='fit')
plt.plot(test2['newevalx'], test2['neweval'], c='b')
plt.plot(test3['newevalx'], test3['neweval'], c='b')
plt.plot(test4['newevalx'], test4['neweval'], c='b')
plt.plot(test5['newevalx'], test5['neweval'], c='b')
plt.plot(test6['newevalx'], test6['neweval'], c='b')
plt.plot(test7['newevalx'], test7['neweval'], c='b')
plt.xlabel('Energy (eV)')
plt.ylabel('Counts per 1 eV bin')
plt.title(str(plottitle))
plt.legend()
plt.show()
plt.close() 


    




