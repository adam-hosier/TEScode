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

theoryres = 10000
minE = 500
maxE = 8000
resStep = (maxE-minE)/theoryres
xray_range = np.arange(minE, maxE, resStep)
EBIT_model_old = mass.materials.filterstack_models['EBIT 2018']
EBIT_model_new = mass.materials.FilterStack(name='2022 Stack')
EBIT_model_new.add_Film(name='Electroplated Au Absorber', material='Au',area_density_g_per_cm2=ufloat(0.00186,0.00006), fill_fraction=ufloat(1.000,0), absorber=True)
EBIT_model_new.add_AlFilmWithOxide(name='50mK Filter', Al_thickness_nm=100)
EBIT_model_new.add_AlFilmWithOxide(name='5K Filter', Al_thickness_nm=112)
EBIT_model_new.add_AlFilmWithOxide(name='50K Filter', Al_thickness_nm=116)
EBIT_model_new.add_Film(name='Ni Mesh', material='Ni', area_density_g_per_cm2=ufloat(0.0134,0.0018), fill_fraction=ufloat(0.170,0.010), absorber=False)
EBIT_model_new.add_LEX_HT(name='Luxel Window #1')

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
ftheoryloc = str('C:\\data\\theory')
teseffloc = str('C:\\data')
theoryEl = 'Nd'
eltitle = theoryEl
#floc = str('C:\\data\\calibratedTES_Dec2022')
ddest = str('C:\\data\\Line_ID_Nd')
date = str('202212')
day = str('21')
runnum = str('0002')

statelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']
statelist = ['A', 'B', 'I', 'O', 'AH']

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
#effdat = pd.read_csv(r""+teseffloc+'\\'+'TES_Efficiency_Dec2022.csv')
#teseff = effdat['Efficiency %']
teseff = dfeff['Efficiency %']
tenergy = tdat['energy']
tintensity = tdat['Intensity']





for s in statelist: 
    state = str(s)
    df[state] = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'.csv')
    #df = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'.csv')

    # counts, bin_edges = np.histogram(df[state]['energy'], bins=numbins, range=(minenergy, maxenergy))
    # df[state+str(' counts')]= counts
    # df[state+str(' bin_edges')] = bin_edges




plt.figure()
plt.plot(df['B']['0'], df['B']['1'])
plt.xlabel('energy eV')
plt.ylabel('photon counts per 1 eV bin')
plt.show()
plt.close() 
asdf
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

state = 'T'
# arry = df[state+str(' counts')]
# arrx = df[state+str(' bin_edges')][:-1]
arry = df[state]['1']
arrx = df[state]['0']
    
    
    
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


res_dict = dict()
a = MultiPeakGaussian(arr = arry, xs = arrx, num_peaks=40, num_poly=3)
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

xfit = res_dict['x0']
yfit = res_dict['y0']
lcycler = cycle(lineSlist)
plt.figure() 
plt.title(str(date)+'_'+str(day)+'_'+str(runnum)+'_'+str(state)+' /// '+str(eltitle)+' '+ str(ebeamen)+' keV , '+str(ebeamcurr)+' mA')
plt.ylabel('Photon counts')
plt.xlabel('Energy (eV)')
plt.minorticks_on() 
plt.plot(arrx, arry, c='b', label='Experimental data')
plt.plot(xfit, yfit, c='g', label='Fitted data')
plt.plot(tx, ty*np.max(arry), c='r', label='Theory (all charge states)')
for k in cstates:
    plt.plot(dfspec['Tx'], dfspec[str(k)+' calcint'],c=np.random.rand(3,),ls=next(lcycler), label='Spectroscopic charge: '+str(k))
plt.legend() 
plt.show() 
plt.close() 