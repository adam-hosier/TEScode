import numpy as np
import csv 
import math as m 
import pandas as pd 
from matplotlib import mlab as Ml 
from scipy.optimize import curve_fit
from scipy import stats
import lmfit 
import scipy.odr.odrpack as odrpack
import matplotlib
import matplotlib.pyplot as plt 
from lmfit import minimize, Parameters, report_fit, Model 
from lmfit.models import GaussianModel
import copy
import os 
from os import listdir 
from os import walk 
from os.path import isfile, join
import glob

##### Need to update this script entirely for TES/x-ray lines for calibration 
floc = str('C:\\Users\\ahosi\\OneDrive\\Desktop\\calibratedTES_Dec2022')
#floc = str('C:\\data\\calibratedTES_Dec2022')
#ddest = str('C:\\data\\Line_ID_Nd')
date = str('202212')
day = str('21')
runnum = str('0002')
calstates = ["A", "I", "O", "AH"]
#statelist = ['T', 'V', 'X', 'Z', 'AB', 'AD', 'AF']
statelist = calstates
coAdd = False
dfall = dict()
minenergy = 500
maxenergy = 5000
binsize = 0.75
numbins = int(np.round((maxenergy-minenergy)/binsize))
findex = np.linspace(0, numbins, num=numbins+1)

for s in statelist: 
    state = str(s)
    dfall[state] = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')
    df = pd.read_csv(r""+floc+'\\'+date+day+'_'+runnum+'_'+state+'photonlist.csv')

    counts, bin_edges = np.histogram(dfall[state]['energy'], bins=numbins, range=(minenergy, maxenergy))
    dfall[state+str(' counts')]= counts
    dfall[state+str(' bin_edges')] = bin_edges

# energy = df['energy']
# time = df['time']

energy = dfall['A']['energy']


counts, bin_edges = np.histogram(energy, bins=numbins, range=(minenergy, maxenergy))
#dattest = np.array((counts, bin_edges), dtype=object).T 
dattest = np.vstack((counts, bin_edges[:-1])).T
df = pd.DataFrame(data=dattest, columns=['counts', 'bin_edges'])

###################################
# plt.figure()
# plt.ylabel('Counts per '+str(binsize)+' eV bin')
# #plt.ylim(bottom=0, top=1.1*np.max(counts))
# #plt.xlim(left=3075, right=3175)
# #plt.xlabel('energy (eV)')
# plt.xlabel('channel number')
# plt.title(date+day+'_'+runnum)
# #plt.title(date+day+'_'+runnum+'_'+state)
# for state in statelist: 
#     #plt.plot(dfall[state+str(' bin_edges')][:-1], dfall[state+str(' counts')], label=state)
#     plt.plot(findex[:-1], dfall[state+str(' counts')], label=state)
# plt.legend()
# plt.show() 
# ##################################


plot = True
#pval = 0.025
pval = 0.16 
cinterval = 95
crit = 100
siglev = 1

fnew = 'C:\\data\\TES_ReCal'

#df = pd.read_csv(r"C:\\Users\\ahosi\\OneDrive\\Desktop\\FileTrans4\\OsIr_NucChargeRadius_ReducedDataTable.csv")

df = df['counts']

#df.drop(df.columns[df.columns.str.contains('Os', case=False)], axis=1, inplace=True)


pixel = findex


# wlen = np.array([
#     1486.2950,
#     1739.3940,
#     2620.8460,
#     3311.1956,
#     4504.9201,
#     6391.0264])


# wlenerr = np.array([0.0100,
#     0.0340,
#     0.0390,
#     0.0060,
#     0.0094,
#     0.0099])


# c = np.array([1314,
#     1652,
#     2828, 
#     3750,
#     5346, 
#     5907])


wlen = np.array([
    1486.2950,
    1739.3940,
    2620.8460,
    3311.1956,
    4504.9201])


wlenerr = np.array([0.0100,
    0.0340,
    0.0390,
    0.0060,
    0.0094])


c = np.array([1314,
    1652,
    2828, 
    3750,
    5346])



r= 10
rb=1 
order=1 
#subtracting last 6 spectra since the CCD plane was moved, requires separate calibration 
num_cols = 1 
m = 0


channel = np.zeros([len(c), num_cols-m], dtype=object)
channel2 = np.zeros([len(c), num_cols-m], dtype=object)
channelerr = np.zeros([len(c), num_cols-m], dtype=object)
color_idx = np.linspace(0,1,num_cols)
resid = np.zeros([len(c), num_cols-m])
sn = np.zeros([len(c), num_cols-m])



def datcoll(df, c, r):  #in single spectra loop
    res = dict()

    y = []   
    ye = []                   

    pix = []

    photc = []                  #photon count


    for i in range(1,2*r):
        y.append(df.iloc[c - r + i])
        pix.append(c-r+i)

        s = y[i-1] 
        photc.append(s)

            
    ye = np.sqrt(photc)

    
    photc = np.array(photc) 
    ye = np.array(ye) 
    pix = np.array(pix) 

    res['y'] = photc 
    res['ye'] = ye 
    res['x'] = pix 
    #res['specnum'] = specn 

    return res 


def fitproc(x, y, ye,c):#in single spectra loop

    fitres = dict()

    def fit(x, A, mu, sig, B):
        return A*np.exp(-(x-mu)**2 / (2 * sig**2))+B

    mod = Model(fit)
    params = Parameters() 
    params.add('A', value = 1000, min=0)
    params.add('mu', value=c)       #, min=c-1, max=c+1
    params.add('sig', value=10, min=0)
    params.add('B', value=1, min=0)

    result = mod.fit(y, params=params, x=x, weights = ye, method='leastsq', nan_policy='omit')
    params.update(result.params)
    xtest = np.linspace(np.min(x), np.max(x), num=1000)
    ytest = mod.eval(x=xtest, params=params)
    # plt.figure() 
    # plt.plot(x, y, label='data')
    # plt.plot(xtest, ytest, label='fit')
    # plt.legend()
    # plt.show()
    # plt.close()
    fitres['A'] = result.params['A'].value
    fitres['Aerr'] = result.params['A'].stderr
    fitres['mu'] = result.params['mu'].value
    fitres['muerr'] = result.params['mu'].stderr
    fitres['sig'] = result.params['sig'].value
    fitres['sigerr'] = result.params['sig'].stderr 
    totN = np.sum(y) 
    fitres['muerr2'] = result.params['sig'].value / np.sqrt(totN)
    fitres['B'] = result.params['mu'].value
    fitres['Berr'] = result.params['mu'].stderr 
    return fitres 


lim = 50

# for k in range(m, num_cols):    

for l in range(len(c)):

    datacollection = datcoll(df, c[l], r) 

    x = datacollection['x']
    y = datacollection['y'] 
    ye = datacollection['ye']
    
    Gaussfit = fitproc(x,y,ye, c[l]) 
    channel[l] = Gaussfit['mu']
    channelerr[l] = Gaussfit['muerr']
    #sn[l,k] = datacollection['specnum']

    resid[l] = c[l] - Gaussfit['mu']

    if Gaussfit['mu'] > c[l] + lim or Gaussfit['mu'] < c[l] - lim:                  #This has to be reconsidered
        channel[l] = 'nan'
        channelerr[l] = 'nan'




#function used to remove NaN's
def nankill(arr, arre, spec):
    
    indexlist = np.zeros([1,0])
    for i in range(len(arr)):
        if arr[i] == 'nan':
            indexlist = np.append(indexlist,i)
        else:
            continue
    
    
    indexlist = np.around(indexlist)
    indexlist = indexlist.astype(int)
    
    arr = np.delete(arr, indexlist)
    arre = np.delete(arre, indexlist) 
    spec = np.delete(spec, indexlist)
    arr = np.array(arr, dtype=np.float64)
    arre = np.array(arre, dtype=np.float64)
    spec = np.array(spec, dtype = np.float64)

    res = dict()
    res['arr'] = arr
    res['arre'] = arre
    res['spec'] = spec 
    return res 

#derivative of the calibration function to calculate confidence bands
def deriv(x, B, C, D, unc):
    
    return (B + 2*C*x + 3*D*x**2)*unc

#partial derivatives w.r.t. fitting coefficients column vector of calibration polynomial to calc confidence bands 
def parts(x):
    
    arr = np.zeros([4,1])

    arr[0,0] = 1
    arr[1,0] = x 
    arr[2,0] = x**2 
    arr[3,0] = x**3


    return arr 

#pixel = np.linspace(0,2047, num=2048)
ap = np.zeros((4,0))

for i in range(len(pixel)):
    ap = np.hstack((ap, parts(pixel[i]) ))



def pixelcal2(cen, cenerr, wlen, wlenerr):#in global loop 
    polyc = dict()
    cb = np.zeros([1,np.shape(ap)[1]])
    cbno = np.zeros([1,np.shape(ap)[1]])
    def pfit(x, A, B, C, D):
        return A + B*x + C*x**2 + D*x**3

    mod = Model(pfit)
    params = Parameters() 
    params.add('A', value=400)
    params.add('B', value=1)
    params.add('C', value=1e-01)
    params.add('D', value=1e-03)

    mod2 = Model(pfit)
    params2 = Parameters() 
    params2.add('A', value=400)
    params2.add('B', value=1)
    params2.add('C', value=1e-01)
    params2.add('D', value=1e-03)


    cen = np.array(cen)
    cenerr = np.array(cenerr)
    wlen = np.array(wlen)
    wlenerr = np.array(wlenerr)
    badcount = 0
    indexlist = np.zeros([1,0])
    for i in range(len(cen)):
        if cen[i] == 'nan':
            indexlist = np.append(indexlist,i)
            badcount += 1
        else:
            continue
    
    
    indexlist = np.around(indexlist)
    indexlist = indexlist.astype(int)
    
    cen = np.delete(cen, indexlist)
    cenerr = np.delete(cenerr, indexlist) 
    wlen = np.delete(wlen, indexlist) 
    wlenerr = np.delete(wlenerr, indexlist) 

    cen = np.array(cen, dtype=np.float64)
    cenerr = np.array(cenerr, dtype=np.float64)
    wlen = np.array(wlen, dtype=np.float64)
    wlenerr = np.array(wlenerr, dtype=np.float64)

    fit = mod.fit(x=cen, data=wlen, params=params, nan_policy='omit', weights = 1/wlenerr, max_nfev=2000)
    params.update(fit.params)
    params2.update(fit.params)
    #print(fit.chisqr)
    def testd(x, B, C, D, chanlu):
        return (B + 2*C*x + 3*D*x**2)*chanlu

    def test(x, A, B, C, D):
        return A + B*x + C*x**2 + D*x**3

    chanluncnm = testd(cen, fit.params['B'].value, fit.params['C'].value, fit.params['D'].value, cenerr)
    
    cutoff = 1.5        #cutoff reduced chi-squared value for iterative systematic uncertainty estimation
    inc = 0.001          #incremental step for " " " " (nm)
    maxit = 1000           #max iterations for " " " "
    b = 0
    pixel = findex
    #print('initial: ', fit.redchi)
    while fit.redchi > cutoff: 

        chanluncnm = testd(cen, fit.params['B'].value, fit.params['C'].value, fit.params['D'].value, cenerr)
        newerr = np.sqrt((wlenerr)**2 + (chanluncnm)**2) + (b*inc)
        fit = mod.fit(x=cen, data=wlen, params=params, nan_policy='omit', weights=1/newerr, max_nfev=2000)
        if b > maxit:
            print('maxit')
            break
        else:
            b+=1

    # if type(sunc) == 'numpy.ndarray':
    #     print('@@@')
    #     sunc = float(sunc)
        
    # else:
    #     sunc = sunc

    #print(type(sunc))
    newerr = np.sqrt((wlenerr)**2 + (chanluncnm)**2) + (b*inc)
    newerrno = np.sqrt((wlenerr)**2 + (chanluncnm)**2)
    fit = mod.fit(x=cen, data=wlen, params=params, nan_policy='omit', weights=1/newerr, max_nfev=2000)
    #print('final: ', fit.redchi)
    params.update(fit.params)
    fitno = mod2.fit(x=cen, data=wlen, params=params2, nan_policy='omit', weights = 1/newerrno, max_nfev=2000)
    params2.update(fitno.params)
    cal = test(pixel, fit.params['A'].value, fit.params['B'].value, fit.params['C'].value, fit.params['D'].value)
    calno = test(pixel, fitno.params['A'].value, fitno.params['B'].value, fitno.params['C'].value, fitno.params['D'].value)
    residual = fit.residual / fit.weights
    residualno = fitno.residual / fitno.weights

    num=len(wlen)
    degf = num - 4
    tval = stats.t.ppf(1-((100-cinterval)/2/100), degf)

    for k in range(np.shape(ap)[1]):

        temp = np.matmul(fit.covar, ap[:,k])
        cb[0,k] = tval*np.sqrt(np.matmul(ap[:,k].T, temp))

        tempno = np.matmul(fitno.covar, ap[:,k])
        cbno[0,k] = tval*np.sqrt(np.matmul(ap[:,k].T, tempno))


    conband = fit.eval_uncertainty(sigma=siglev, x=cal)
    conband = np.reshape(conband, (np.shape(conband)[0],1)).T
    conbandno = fitno.eval_uncertainty(sigma=siglev, x=cal)

  
    polyc['cbno'] = cbno 
    polyc['paramsno'] = fitno.params
    polyc['covarno'] = fitno.covar
    polyc['residualno'] = residualno
    polyc['newerrno'] = newerrno
    polyc['cb2'] = conband
    polyc['cb2no'] = conbandno
    polyc['cb'] = cb
    polyc['covar'] = fit.covar
    polyc['K0'] = fit.params['A'].value
    polyc['K0e'] = fit.params['A'].stderr
    polyc['K1'] = fit.params['B'].value
    polyc['K1e'] = fit.params['B'].stderr
    polyc['K2'] = fit.params['C'].value
    polyc['K2e'] = fit.params['C'].stderr
    polyc['K3'] = fit.params['D'].value
    polyc['K3e'] = fit.params['D'].stderr
    polyc['red_chi_sq'] = fit.redchi
    polyc['badcount'] = badcount
    polyc['calibration'] = cal
    polyc['sysunc'] = b*inc
    polyc['newerr'] = newerr
    polyc['residual'] = residual
    polyc['tval'] = tval

    return polyc


polyco = np.zeros([4,num_cols], dtype=object)
polycoerr = np.zeros([4,num_cols], dtype=object)
badc = np.zeros([1,num_cols])


#specnumber = sn[0,:]

#polynomial for temporal evolution of lines 
def fun(x, A, B, C, D, E, F):
    return A + B*x + C*x**2 + D*x**3 + E*x**4 + F*x**5

#histogram analysis to determine the distribution of the data around the fit 
def Histo(resid, bins, bin_range):
    binneddata = np.histogram(resid, bins=bins, range=[-bin_range, bin_range])         #binned data 


    HisRes = dict()
    yuncert = np.zeros([1, len(binneddata[0][:])])

    uncerterr = np.zeros([1, len(binneddata[0][:])])

    freq = np.zeros([1, len(binneddata[0][:])])
 
    totalobs = 0


    for i in range(len(binneddata[1][:])-1):
        
        yuncert[0,i] = (1/2) * (binneddata[1][i] + binneddata[1][i+1])

        if binneddata[0][i] == 0:
            uncerterr[0,i] = 'nan'
        else:
            uncerterr[0,i] = 1 / np.sqrt(binneddata[0][i])

        
        freq[0,i] = (binneddata[0][i])

        totalobs += binneddata[0][i] 


    HisRes['totalobs'] = totalobs 


    mod1 = GaussianModel()

    params = mod1.make_params(amp = 10, cen = 0, sig = 0.1)

    errfit = mod1.fit(freq, params, x=yuncert, nan_policy='omit', weights = uncerterr)  #weights = uncerterr

    params.update(errfit.params)


    plt.figure()

    for i in range(bins):
        
        plt.scatter(yuncert[0,i], freq[0,i], color = 'b')

    xgraph = np.linspace(-bin_range, bin_range, num=1000)

    values1 = params.valuesdict()


    cen1 = values1['center']

    ###################################Graphing
    
    resideval = errfit.eval(params, x=xgraph)

    plt.plot(xgraph, resideval, color='c')

    plt.ylabel('Frequency (Counts)')
    plt.xlabel('Binned Residuals (Unweighted, eV)')
    plt.axvline(x=cen1, color='b')

    plt.close()
    #plt.show()
    
    ####################################
    calcparams = params.valuesdict()

    HisRes['sigma'] = calcparams['sigma']


    return HisRes


def OutlierDet(cen, cenerr, wlen, wlenerr, crit):       #single point removal
    #badlines = np.zeros([1,0])
    sunc2 = np.zeros([1,0])
    newres = np.zeros([1,0])
    for i in range(len(cen)):
        arr = copy.deepcopy(cen)            #making copy of arrays to avoid altering the original arrays
        arrerr = copy.deepcopy(cenerr)
        warr = copy.deepcopy(wlen)
        warrerr = copy.deepcopy(wlenerr)

        arr2 = np.delete(arr, i)             #deleting ith row in the data (rotating through all data for single removal)
        arrerr2 = np.delete(arrerr, i)
        warr2 = np.delete(warr, i)
        warrerr2 = np.delete(warrerr, i)

        newfit = pixelcal2(arr2, arrerr2, warr2, warrerr2)
        sunc2 = np.append(sunc2, newfit['sysunc'])
        val = newfit['K0'] + (arr[i])*newfit['K1'] + (newfit['K2']) * (arr[i])**2 + (newfit['K3']) * (arr[i])**3 
        newres = np.append(newres, np.abs(warr[i]-val))

    maxres = np.max(newres)
    maxin = np.argmax(newres)
    #print(np.around(arr[maxin]).astype(int))
    #print(newfit['confband'][0,np.around(arr[maxin]).astype(int)])
    if (maxres/(newfit['cb'][0,np.around(arr[maxin]).astype(int)])) >= crit:            #rejection criteria based on 
        #badlines = np.append(badlines, warr[maxin]) 

        badlines = warr[maxin]                                         #only the distance from the confidenceband
        newarr = np.delete(cen, maxin)
        newarrerr = np.delete(cenerr, maxin)
        newwarr = np.delete(wlen, maxin)
        newwarrerr = np.delete(wlenerr, maxin)
    
    if (maxres/(newfit['cb'][0,np.around(arr[maxin]).astype(int)])) < crit:
        #badlines = np.append(badlines, None)

        badlines = None
        newarr = cen
        newarrerr = cenerr
        newwarr = wlen
        newwarrerr = wlenerr
        newfit = pixelcal2(cen, cenerr, wlen, wlenerr)

    else: 
        #somehow going to here
        badlines = None
        newarr = cen
        newarrerr = cenerr
        newwarr = wlen
        newwarrerr = wlenerr
        newfit = pixelcal2(cen, cenerr, wlen, wlenerr)

    pixel = findex
    ODcal = newfit['K0'] + newfit['K1']*pixel + newfit['K2']*pixel**2 + newfit['K3']*pixel**3

    ODres = {}
    ODres['arr'] = newarr
    ODres['arrerr'] = newarrerr
    ODres['wlen'] = newwarr
    ODres['wlenerr'] = newwarrerr 
    ODres['badlines'] = badlines 
    ODres['sysunc'] = newfit['sysunc']
    ODres['K0'] = newfit['K0']
    ODres['K1'] = newfit['K1']
    ODres['K2'] = newfit['K2'] 
    ODres['K3'] = newfit['K3'] 
    ODres['K0e'] = newfit['K0e']
    ODres['K1e'] = newfit['K1e']
    ODres['K2e'] = newfit['K2e']
    ODres['K3e'] = newfit['K3e']
    ODres['calibration'] = ODcal
    ODres['residual'] = newfit['residual']
    ODres['newerr'] = newfit['newerr']
    ODres['confband'] = newfit['cb']
    
    return ODres
    

def nankill2(arr):
    rows, columns = np.shape(arr)
    print('@@@@@@@@@@@@')
    i = 0
    a=0
    for i in range(rows):
        if np.any(arr[i,:]==np.nan) or np.any(arr[i,:]=='nan') or np.any(arr[i,:]==None):
            arr = np.delete(arr, i, axis=0)
            a+=1
            break
        else:
            continue
    
    arr = np.reshape(arr, (rows-a, columns))
    return arr



bins=12

Kpars = ['K0', 'K1', 'K2', 'K3']

####Actual calibration procedure/fitting 
confband = np.zeros([num_cols, np.shape(ap)[1]])
names = []
calib = np.zeros([num_cols, len(pixel)])
sysunc = np.zeros([num_cols,1])
sysunc2 = np.zeros([num_cols,1])
calnum = []
cali = []
tval1 = []
for i in range(num_cols):

    names.append(str('test'))
    data1 = channel[:,i]    #Ne data

    data = data1
    err1 = channelerr[:,i]

    err = err1
    wavelerr = wlenerr
    wavel = wlen
    arr = np.array((data, err, wavel, wavelerr)).T


    while np.any(arr=='nan'):
        print('flag')
        arr = nankill2(arr)

    arr = arr[arr[:,0].argsort()]
    pcal = pixelcal2(arr[:,0], arr[:,1], arr[:,2], arr[:,3])

    # plt.figure() 
    # plt.scatter(arr[:,2], arr[:,0])
    # plt.show() 
    # plt.close() 

    polyco[0,i] = pcal['K0']
    polyco[1,i] = pcal['K1']
    polyco[2,i] = pcal['K2']
    polyco[3,i] = pcal['K3']
    polycoerr[0,i] = pcal['K0e']
    polycoerr[1,i] = pcal['K1e']
    polycoerr[2,i] = pcal['K2e']
    polycoerr[3,i] = pcal['K3e']
    sysunc[i,0] = pcal['sysunc']
    #confband[i,:] = pcal['confband']
    confband[i,:] = pcal['cb2']
    calib[i,:] = pcal['calibration']
    
    # print(pcal['cb2'])
    # print(pcal['confband'])
    # print('@@@@')
    cal = calib[i,:]
    cband = confband[i,:]
    preparams = np.array([pcal['K0'], pcal['K1'], pcal['K2'], pcal['K3']])
    preparerr = np.array([pcal['K0e'], pcal['K1e'], pcal['K2e'], pcal['K3e']])
    prepar = np.array([Kpars, preparams, preparerr]).T
    prepardf = pd.DataFrame(data=prepar, columns=['Parameters', 'Value', 'Uncertainty'])
    # prepardf.to_csv(fnew  + '/'+ 'Cal'+str(i+1) + names[i]+'_' + 'PolyParameters_Pre_Removal'+ '.csv', index=False)

    precal = np.reshape(cal, (np.shape(cal)[0],))
    preconfband = np.reshape(cband, (np.shape(cband)[0],))
    preres = np.array([precal, preconfband]).T
    preresdf = pd.DataFrame(data=preres, columns=['calibration', 'conf band'])
    # preresdf.to_csv(fnew  + '/'+ 'Cal'+str(i+1) + names[i]+'_' + 'Calibration_Pre_Removal'+ '.csv', index=False)



    totbad = []
    a = arr[:,0]
    b = arr[:,1]
    c = arr[:,2] 
    d = arr[:,3]


    count = 0
    OutDet = OutlierDet(a, b, c, d,crit)
    print('Beginning outlier detection...')
    while OutDet['badlines'] is not None:
        
        a2 = OutDet['arr']
        b2 = OutDet['arrerr']
        c2 = OutDet['wlen']
        d2 = OutDet['wlenerr']
        totbad.append(OutDet['badlines'])
        sysunc2[i,0] = OutDet['sysunc']

        a = a2 
        b = b2 
        c = c2 
        d = d2 
        
        if count > np.shape(arr)[0]:
            print('max looped')
            break
        else:
            count +=1 
            OutDet = OutlierDet(a, b, c, d,crit)
    print('Finished outlier detection')
    print(totbad)
    #OutDet = OutlierDet(a, b, c, d,crit)
    newfit = pixelcal2(OutDet['arr'], OutDet['arrerr'], OutDet['wlen'], OutDet['wlenerr'])
    postparams = np.array([newfit['K0'], newfit['K1'], newfit['K2'], newfit['K3']])
    postparerr = np.array([newfit['K0e'], newfit['K1e'], newfit['K2e'], newfit['K3e']])
    postpar = np.array([Kpars, postparams, postparerr]).T
    postpardf = pd.DataFrame(data=postpar, columns=['Parameters', 'Value', 'Uncertainty'])
    postpardf.to_csv(fnew  +  names[i]+'_' + 'PolyParameters'+ '.csv', index=False)
    
    totbad = np.reshape(totbad, (len(totbad),1))
    dfz2 = totbad 
    newdfz2 = pd.DataFrame(data=dfz2, columns=['lines removed (nm)'])

    newdfz2.to_csv(fnew  + '/' +'Cal'+str(i+1) + names[i]+'_' + 'RemovedLines'+ '.csv', index=False)

    sysunc2[i,0] = newfit['sysunc']
    tval1.append(newfit['tval'])
    if plot==True:
        # plt.figure() 
        # plt.title(str(df.columns[i]) + '  ,  ' + str(dfb.columns[i]))
        # plt.xlabel('Wavelength (nm)')
        # plt.ylabel('Residual (unweighted)')
        # plt.scatter(x=arr[:,2], y=pcal['residual'], c='b', marker="^", label='Residual')
        # plt.scatter(x=c, y=newfit['residual'], c='r', marker="^", label='Residual(removed outliers)')
        # plt.errorbar(x=arr[:,2], y=pcal['residual'], yerr=pcal['newerr'], ecolor='b', ls='none')
        # plt.errorbar(x=c, y=newfit['residual'], yerr=newfit['newerr'], ecolor='r', ls='none')
        # plt.ylim([-10**(-2), 10**(-2)])
        # plt.xlim([np.min(arr[:,2]), np.max(arr[:,2])])
        # plt.ylim([-4*10**(-3), 4*10**(-3)])
        # plt.plot(pcal['calibration'], confband[i,:], c='b', label='95% conf band')
        # plt.plot(pcal['calibration'], -confband[i,:], c='b')
        # plt.plot(newfit['calibration'], newfit['confband'][0,:], c='r', label='95% conf band(removed outliers)')
        # plt.plot(newfit['calibration'], - newfit['confband'][0,:], c='r')
        # plt.legend()
        # plt.savefig(fnew + '/'  + 'Cal'+str(i+1)+ names[i]+ '.pdf')

        
        # #plt.show()
        
        # plt.close()

        plt.figure(figsize=(15,9)) 
        plt.title(str('test'))
        plt.xlabel('Energy (eV)')
        plt.ylabel('Residual (unweighted)')
        #plt.scatter(x=arr[:,2], y=pcal['residual'], c='b', marker="^", label='Residual')
        plt.scatter(x=c, y=newfit['residual'], c='r', marker="^", label='Residual(with systematic unc: '+str(newfit['sysunc'])+')')
        plt.scatter(x=c, y=newfit['residualno'], c='b', marker="^", label='Residual (without)')
        #plt.errorbar(x=arr[:,2], y=pcal['residual'], yerr=pcal['newerr'], ecolor='b', ls='none')
        plt.errorbar(x=c, y=newfit['residual'], yerr=newfit['newerr'], ecolor='r', ls='none')
        plt.errorbar(x=c, y=newfit['residualno'], yerr=newfit['newerrno'], ecolor='b', ls='none')
        #plt.ylim([-10**(-2), 10**(-2)])
        #plt.xlim([np.min(arr[:,2]), np.max(arr[:,2])])
        #plt.xlim([4.00, 20.00])
        #plt.ylim([-5*10**(-3), 5*10**(-3)])
        #plt.plot(pcal['calibration'], confband[i,:], c='b', label='95% conf band')
        #plt.plot(pcal['calibration'], -confband[i,:], c='b')
        #plt.plot(newfit['calibration'], newfit['cb2'][0,:], c='r', label=str(siglev)+'-sigma conf band(with)')
        #plt.plot(newfit['calibration'], - newfit['cb2'][0,:], c='r')
        #plt.plot(newfit['calibration'], newfit['cb2no'][0,:], c='b', label=str(siglev)+'-sigma conf band (without)')
        plt.plot(newfit['calibration'], newfit['cb'][0,:], c='r', label=str(cinterval)+'% conf band(with)')
        plt.plot(newfit['calibration'], - newfit['cb'][0,:], c='r')
        plt.plot(newfit['calibration'], newfit['cbno'][0,:], c='b', label=str(cinterval)+'% conf band (without)')
        plt.minorticks_on()
        plt.plot(newfit['calibration'], -newfit['cbno'][0,:], c='b')
        plt.legend()
        plt.savefig(fnew + names[i]+ '.pdf')
        #plt.show()
        
        plt.close()


        ###spectra of Ne and Bg together 
        plt.figure(figsize=(22,9)) 
        plt.title(str('test'))
        plt.xlabel('Energy (eV)')
        plt.ylabel('ADU (arb)')
        plt.plot(newfit['calibration'][:-1], df, c='r', label=df.keys()[i])
        #plt.plot(newfit['calibration'], dfb[dfb.keys()[i]], c='b', label=dfb.keys()[i]) 
        plt.xlim((np.min(newfit['calibration']),np.max(newfit['calibration'])) )
        for l in range(len(wlen)):
            plt.axvline(x=wlen[l], ls='--', c='tab:brown')

        plt.legend() 
        plt.minorticks_on()
        plt.savefig(fnew+'_Spectra_'+ names[i]+ '.pdf')

        plt.close()
        

    
    cband = np.reshape(newfit['cb'], (np.shape(newfit['cb'])[1],))
    dfz = np.array([newfit['calibration'], cband]).T

    newdfz = pd.DataFrame(data=dfz, columns=['calibration', 'confband'])
    #newdfz.to_csv(fnew  + '/' +'Cal'+str(i+1)+ names[i]+'_' + 'Calibration'+ '.csv', index=False)
    newdfz.to_csv(fnew  +  names[i]+'_' + 'Calibration'+ '.csv', index=False)
    calnum.append('Cal'+str(i+1))
    cali.append(newfit['calibration'])

    OsNa = parts(567)
    OsMg = parts(597)
    IrNa = parts(544)
    IrMg = parts(573)
    '''
    temp1 = np.matmul(newfit['covar'], OsNa[:,0])
    temp = 2*np.sqrt(np.matmul(OsNa[:,0].T, temp1))
    print('OsNa: ', temp) 

    temp1 = np.matmul(newfit['covar'], OsMg[:,0])
    temp = 2*np.sqrt(np.matmul(OsMg[:,0].T, temp1))
    print('OsMg: ', temp) 

    temp1 = np.matmul(newfit['covar'], IrNa[:,0])
    temp = 2*np.sqrt(np.matmul(IrNa[:,0].T, temp1))
    print('IrNa: ', temp) 

    temp1 = np.matmul(newfit['covar'], IrMg[:,0])
    temp = 2*np.sqrt(np.matmul(IrMg[:,0].T, temp1))
    print('IrMg: ', temp) 
    '''
    dfcovar = np.array(newfit['covar'])
    #print(newfit['covar']) 
    #print(dfcovar)
    newdfcovar = pd.DataFrame(data=dfcovar, columns=['K0', 'K1', 'K2', 'K3'])
    #newdfcovar.to_csv(fnew  + '/' + 'Cal'+str(i+1)+ names[i]+'_' + 'CovarianceMatrix'+ '.csv', index=False)
    newdfcovar.to_csv(fnew + names[i]+'_' + 'CovarianceMatrix'+ '.csv', index=False)
    print('Calibration #'+ str(i+1)+' completed.')


sysunc = np.reshape(sysunc, (np.shape(sysunc)[0],))
sysunc2 = np.reshape(sysunc2, (np.shape(sysunc2)[0],))
tval1 = np.reshape(tval1, (np.shape(tval1)[0],))
dfz3 = np.array([sysunc, sysunc2, tval1]).T
newdfz3 = pd.DataFrame(data=dfz3, columns=['Before removal', 'after removal', 'tval'])
newdfz3.to_csv(fnew + 'Systematic_Uncertainty.csv', index=False)

cali = np.array(cali).T
newcali = pd.DataFrame(data=cali, columns=calnum)
newcali.to_csv(fnew+'TEStestCal.csv', index=False)