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

plot = True
#pval = 0.025
pval = 0.16 
cinterval = 95
crit = 3.5
siglev = 1

fnew = 'C:\\Users\\ahosi\\Desktop\\Calibration6b'

df = pd.read_csv(r"C:\\Users\\ahosi\\OneDrive\\Desktop\\FileTrans4\\OsIr_NucChargeRadius_ReducedDataTable.csv")


df.drop(df.columns[df.columns.str.contains('Os', case=False)], axis=1, inplace=True)


pixel = np.linspace(1,2048, num=2048)


wlen = np.array([
#    6.73832,
#    6.894,
    7.5764,
    8.80929,
    9.7502,
#    9.81166,
#    9.82629,
#    10.2909,
#    10.3084,
##    10.6138,
#    11.1136,
#    11.42,
    11.6691,
#    12.2516,
#    12.2702,
    12.7676,
    14.3314,
    14.7138,
    19.5004])
#    19.6233,
#    19.6526])

wlenerr = np.array([
#    0.00012,
#    0.002,
    0.0004,
    0.00014,
    0.0004,
#    0.00021,
#    0.00021,
#    0.00019,
#    0.00019,
##    0.0005,
#    0.0018,
#    0.004,
    0.0005,
#    0.00017,
#    0.00017,
    0.0007,
    0.0007,
    0.0007,
    0.0008])
#    0.00042,
#    0.00042])

c = np.array([
#    459,
#    482,
    586,
    762,
    891,
#    899,
#    901,
#    962,
#    964,
##    1004,
#    1068,
#    1106,
    1137,
#    1208,
#    1210,
    1270,
    1451,
    1494,
    1999])
#    2011,
#    2014])




r= 4
rb=1 
order=1 
#subtracting last 6 spectra since the CCD plane was moved, requires separate calibration 
num_cols = len(df.columns) 
m = 0


channel = np.zeros([len(c), num_cols-m], dtype=object)
channel2 = np.zeros([len(c), num_cols-m], dtype=object)
channelerr = np.zeros([len(c), num_cols-m], dtype=object)
color_idx = np.linspace(0,1,num_cols)
resid = np.zeros([len(c), num_cols-m])
sn = np.zeros([len(c), num_cols-m])



def datcoll(df, df2, c, r, rb, order, num_cols, k):  #in single spectra loop
    res = dict()

    y = []   
    ye = []                   

    pix = []
    bg = 0                  #background outside of radius of data 
    photc = []                  #photon count

    specn = df2.iloc[0,k]
    for i in range(rb):
        #bg += (1/(2*rb))*(df.iloc[c-r-i,k] + df.iloc[c+r+i,k])
        bg = np.min(df.iloc[:,k])

    #print(np.min(df.iloc[:,k]))
    #bg = 1
    #Collecting data around a point of interest defined by 'c' and 'r' above, along with centroid calcs
    for i in range(1,2*r):
        y.append(df.iloc[c - r + i, k])
        
        pix.append(c-r+i)
        #if y[i-1] > bg:
        s = y[i-1]-bg    
        photc.append(s)
            
        #else:                           #if background is higher than selected data point 
            ##photc.append(s)
            
    ye = np.sqrt(photc)

    
    photc = np.array(photc) 
    ye = np.array(ye) 
    pix = np.array(pix) 
    #plt.figure()
    #plt.plot(pix,photc, color=plt.cm.cool(color_idx[k]))

    #plt.ylabel(' intensity (ADU)')
    #plt.xlabel('channel number')
    res['y'] = photc 
    res['ye'] = ye 
    res['x'] = pix 
    res['specnum'] = specn 
    #plt.show()
    return res 


def fitproc(x, y, ye,c):#in single spectra loop

    fitres = dict()

    def fit(x, A, mu, sig):
        return A*np.exp(-(x-mu)**2 / (2 * sig**2)) 

    mod = Model(fit)
    params = Parameters() 
    params.add('A', value = 50, min=0)
    params.add('mu', value=c)       #, min=c-1, max=c+1
    params.add('sig', value=1, min=0)
    

    result = mod.fit(y, params=params, x=x, weights = None, method='leastsq', nan_policy='omit')
    params.update(result.params)

    fitres['A'] = result.params['A'].value
    fitres['Aerr'] = result.params['A'].stderr
    fitres['mu'] = result.params['mu'].value
    fitres['muerr'] = result.params['mu'].stderr
    fitres['sig'] = result.params['sig'].value
    fitres['sigerr'] = result.params['sig'].stderr 
    totN = np.sum(y) 
    fitres['muerr2'] = result.params['sig'].value / np.sqrt(totN)

    return fitres 


lim = 50
###Data collection across the Ne and Bg spectra (correct)
for k in range(m, num_cols):    

    for l in range(len(c)):

        datacollection = datcoll(df, df2, c[l], r, rb, order, num_cols, k) 
 
        x = datacollection['x']
        y = datacollection['y'] 
        ye = datacollection['ye']
        
        Gaussfit = fitproc(x,y,ye, c[l]) 
        channel[l,k] = Gaussfit['mu']
        channelerr[l,k] = Gaussfit['muerr']
        sn[l,k] = datacollection['specnum']

        resid[l,k] = c[l] - Gaussfit['mu']

        if Gaussfit['mu'] > c[l] + lim or Gaussfit['mu'] < c[l] - lim:                  #This has to be reconsidered
            channel[l,k] = 'nan'
            channelerr[l,k] = 'nan'
            #channel[l,k] = np.nan
            #channelerr[l,k] = np.nan



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

pixel = np.linspace(0,2047, num=2048)
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
    params.add('A', value=4)
    params.add('B', value=1e-02)
    params.add('C', value=1e-04)
    params.add('D', value=1e-07)

    mod2 = Model(pfit)
    params2 = Parameters() 
    params2.add('A', value=4)
    params2.add('B', value=1e-02)
    params2.add('C', value=1e-04)
    params2.add('D', value=1e-07)


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
    inc = 0.000005          #incremental step for " " " " (nm)
    maxit = 1000           #max iterations for " " " "
    b = 0
    pixel = np.linspace(0, 2047, num=2048)
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


specnumber = sn[0,:]

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


    pixel = np.linspace(0, 2047, num=2048)
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

    names.append(str(df.columns[i])+str(dfb.columns[i]))
    data1 = channel[:,i]    #Ne data

    data = np.append(data1, data2)
    err1 = channelerr[:,i]

    err = np.append(err1, err2)
    wavelerr = np.append(wlenerr, wlenerrbg)
    wavel = np.append(wlen, wlenbg)
    arr = np.array((data, err, wavel, wavelerr)).T


    while np.any(arr=='nan'):
        print('flag')
        arr = nankill2(arr)

    arr = arr[arr[:,0].argsort()]
    pcal = pixelcal2(arr[:,0], arr[:,1], arr[:,2], arr[:,3])
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
    prepardf.to_csv(fnew  + '/'+ 'Cal'+str(i+1) + names[i]+'_' + 'PolyParameters_Pre_Removal'+ '.csv', index=False)

    precal = np.reshape(cal, (np.shape(cal)[0],))
    preconfband = np.reshape(cband, (np.shape(cband)[0],))
    preres = np.array([precal, preconfband]).T
    preresdf = pd.DataFrame(data=preres, columns=['calibration', 'conf band'])
    preresdf.to_csv(fnew  + '/'+ 'Cal'+str(i+1) + names[i]+'_' + 'Calibration_Pre_Removal'+ '.csv', index=False)



    totbad = []
    a = arr[:,0]
    b = arr[:,1]
    c = arr[:,2] 
    d = arr[:,3]


    count = 0
    OutDet = OutlierDet(a, b, c, d,crit)

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
    print(totbad)
    #OutDet = OutlierDet(a, b, c, d,crit)
    newfit = pixelcal2(OutDet['arr'], OutDet['arrerr'], OutDet['wlen'], OutDet['wlenerr'])
    postparams = np.array([newfit['K0'], newfit['K1'], newfit['K2'], newfit['K3']])
    postparerr = np.array([newfit['K0e'], newfit['K1e'], newfit['K2e'], newfit['K3e']])
    postpar = np.array([Kpars, postparams, postparerr]).T
    postpardf = pd.DataFrame(data=postpar, columns=['Parameters', 'Value', 'Uncertainty'])
    postpardf.to_csv(fnew  + '/' +'Cal'+str(i+1)+ '/' + names[i]+'_' + 'PolyParameters'+ '.csv', index=False)
    
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
        plt.title(str(df.columns[i]) + '  ,  ' + str(dfb.columns[i]))
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Residual (unweighted)')
        #plt.scatter(x=arr[:,2], y=pcal['residual'], c='b', marker="^", label='Residual')
        plt.scatter(x=c, y=newfit['residual'], c='r', marker="^", label='Residual(with systematic unc: '+str(newfit['sysunc'])+')')
        plt.scatter(x=c, y=newfit['residualno'], c='b', marker="^", label='Residual (without)')
        #plt.errorbar(x=arr[:,2], y=pcal['residual'], yerr=pcal['newerr'], ecolor='b', ls='none')
        plt.errorbar(x=c, y=newfit['residual'], yerr=newfit['newerr'], ecolor='r', ls='none')
        plt.errorbar(x=c, y=newfit['residualno'], yerr=newfit['newerrno'], ecolor='b', ls='none')
        #plt.ylim([-10**(-2), 10**(-2)])
        #plt.xlim([np.min(arr[:,2]), np.max(arr[:,2])])
        plt.xlim([4.00, 20.00])
        plt.ylim([-5*10**(-3), 5*10**(-3)])
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
        plt.savefig(fnew2 + '/'  + 'Cal'+str(i+1)+'_'+ names[i]+ '.pdf')
        #plt.show()
        
        plt.close()


        ###spectra of Ne and Bg together 
        plt.figure(figsize=(22,9)) 
        plt.title(str(df.columns[i]) + '  ,  ' + str(dfb.columns[i]))
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('ADU (arb)')
        plt.plot(newfit['calibration'], df[df.keys()[i]], c='r', label=df.keys()[i])
        plt.plot(newfit['calibration'], dfb[dfb.keys()[i]], c='b', label=dfb.keys()[i]) 
        plt.xlim((np.min(newfit['calibration']),np.max(newfit['calibration'])) )
        for l in range(len(wlen)):
            plt.axvline(x=wlen[l], ls='--', c='tab:brown')

        for l in range(len(wlenbg)):
            plt.axvline(x=wlenbg[l], ls='--', c='tab:gray')

        plt.legend() 
        plt.minorticks_on()
        plt.savefig(fnew2+'/'+'Cal'+str(i+1)+'_Spectra_'+ names[i]+ '.pdf')

        plt.close()
        

    
    cband = np.reshape(newfit['cb'], (np.shape(newfit['cb'])[1],))
    dfz = np.array([newfit['calibration'], cband]).T

    newdfz = pd.DataFrame(data=dfz, columns=['calibration', 'confband'])
    #newdfz.to_csv(fnew  + '/' +'Cal'+str(i+1)+ names[i]+'_' + 'Calibration'+ '.csv', index=False)
    newdfz.to_csv(fnew  + '/' +'Cal'+str(i+1)+ '/'+ names[i]+'_' + 'Calibration'+ '.csv', index=False)
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
    newdfcovar.to_csv(fnew  + '/' + 'Cal'+str(i+1)+ '/'+ names[i]+'_' + 'CovarianceMatrix'+ '.csv', index=False)
    print('Calibration #'+ str(i+1)+' completed.')

sysunc = np.reshape(sysunc, (np.shape(sysunc)[0],))
sysunc2 = np.reshape(sysunc2, (np.shape(sysunc2)[0],))
tval1 = np.reshape(tval1, (np.shape(tval1)[0],))
dfz3 = np.array([sysunc, sysunc2, tval1]).T
newdfz3 = pd.DataFrame(data=dfz3, columns=['Before removal', 'after removal', 'tval'])
newdfz3.to_csv(fnew + 'Systematic_Uncertainty.csv', index=False)

cali = np.array(cali).T
newcali = pd.DataFrame(data=cali, columns=calnum)
newcali.to_csv(fnew+'OsIrWavelengthCal6b.csv', index=False)