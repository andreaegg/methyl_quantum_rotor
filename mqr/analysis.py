from cvxopt import normal
from matplotlib.cbook import ls_mapper
import matplotlib.pyplot as plt
import numpy as np
import scipy as spy
import math as m
import deerlab as dl

def baselinecorr(B,int,points=100,side='right',degree=1,mode='subtraction'):

    int = np.real(int)

    if side is 'left':
        param = np.polyfit(B[0:points+1],int[0:points+1],degree)
        baseline = np.polyval(param,B)
    elif side is 'right':
        param = np.polyfit(B[len(B)-points:len(B)+1],int[len(int)-points:len(int)+1],degree)
        baseline = np.polyval(param,B)
    elif side is 'both':
        param = np.polyfit([B[0:points+1],B[len(B)-points:len(B)+1]],[int[0:points+1],int[len(int)-points:len(int)+1]],degree)
        baseline = np.polyval(param,B)

    if mode is 'division':
        corrint = int/baseline
    elif mode is 'subtraction':
        corrint = int -baseline

    return

def phasecorr(trace):

    def rms_phasecorr(trace,phi):
        imagtr = np.imag(trace*np.exp(-1j*phi))
        rms    = np.sqrt(1/len(trace))*np.sum(imagtr*imagtr)
        return rms

    phi0 = np.arctan2(np.imag(trace[len(trace)-1]),np.real(trace[len(trace)-1]))
    model = lambda x: rms_phasecorr(trace,x)
    phi  = spy.optimize.fmin(model,phi0,disp=False)
    trcorr = trace*np.exp(-1j*phi)
    return trcorr

def stitch_ESEEMdataset(ts,sigs,phase=True,offset=False,normalization=True,plotflag=False):

    if phase is True:
        for k in range(len(sigs)):
            sigs[k] = phasecorr(sigs[k])

    if type(ts) is list:
        if type(sigs) is list:
            t    = ts[0]
            data = sigs[0]
            for k in range(1,len(ts)):
                # find overlap
                idx = np.nonzero(ts[k] == ts[k-1][len(ts[k-1])-1])[0]
                if len(idx) == 0:
                    idx = np.argmin(abs(ts[k-1][len(ts[k-1])-1] - ts[k]))
                else:
                    idx = idx[0]# offset correction
                if offset is True:
                    diff = np.mean(np.real(sigs[k-1][len(sigs[k-1])-idx-1:len(sigs[k-1])])) - np.mean(np.real(sigs[k][0:idx+1]))
                    sigs[k] = sigs[k] + diff
                t   = np.concatenate((t,ts[k][idx+1:len(ts[k])]))
                data = np.concatenate((data,sigs[k][idx+1:len(sigs[k])]))


    # optional normalization
    if normalization is True:
        data = data/np.amax(data)      

        if plotflag is True:
            for k in range(len(ts)):
                if normalization is True:
                    plt.plot(ts[k],np.real(sigs[k])/np.real(sigs[0][0]),'k')
                    plt.plot(ts[k],np.imag(sigs[k])/np.real(sigs[0][0]),'k')
                else:
                    plt.plot(ts[k],np.real(sigs[k]),'k')
                    plt.plot(ts[k],np.imag(sigs[k]),'k')
                
            plt.plot(t,np.real(data),'b')
            plt.plot(t,np.imag(data),'b')
            plt.xlim(0,np.amax(t))
            plt.xlabel('T [$\mu$s]')
            plt.ylabel('intensity [a.u.]')
            plt.legend('raw','stitched,corrected')

    return t,data


def adjust_ussignal_kernel(t,tK,signal,K):

    factor = np.round((tK[1]-tK[0])/(t[1]-t[0]),decimals=1)
    lK = K.shape[1]
    lS = len(signal)

    if lK < lS:
        currt = t[0:lK]
        if (t[0] < tK[0]):
            for n in range(lK):
                start = np.nonzero(np.round(currt,decimals=3) == np.round(tK[n],decimals=3))
                if start:
                    break
            start  = start[0][0]
            t      = t[start:lS]
            signal = signal[start:lS]
            tK     = tK[n:lK]
            K      = K[:,n:lK]
        elif t[0] > tK[0]:
            for n in range(lK):
                start = np.nonzero(np.round(tK,decimals=3) == np.round(currt[n],decimals=3))
                if start:
                    break
            start  = start[0][0]
            tK     = tK[start:lK]
            K      = K[:,start:lK]
            t      = tK[n:lS]
            signal = signal[n:lS]
    
    elif lK > lS:
        currK = tK[0:lS]
        if t[0] < tK[0]:
            for n in range(lS):
                start = np.nonzero(np.round(t,decimals=3) == np.round(currK[n],decimals=3))
                if start:
                    break
            start  = start[0][0]
            t      = t[start:lS]
            signal = signal[start:lS]
            tK     = tK[n:lK]
            K      = K[:,n:lK]
        elif t[0] > tK[0]:
            for n in range(lS):
                start = np.nonzero(np.round(currK,decimals=3) == np.round(t[n],decimals=3))
                if start:
                    break
            start  = start[0][0]
            tK     = tK[start:lK]
            K      = K[:,start:lK]
            t      = t[n:lS]
            signal = signal[n:lS]
    else:
        if t[0] < tK[0]:
            start = np.nonzero(np.round(t,decimals=3) == np.round(tK,decimals=3))
            t      = t[start:lS]
            signal = signal[start:lS]
            tK     = tK[start:lK]
            K      = K[:,start:lK]
        elif t[0] > tK[0]:
            start = np.nonzero(np.round(tK,decimals=3) == np.round(t,decimals=3))
            tK     = tK[start:lK]
            K      = K[:,start:lK]
            t      = tK[start:lS]
            signal = signal[start:lS]

    if factor == 1:
        if len(t) > len(tK):
            ende = np.nonzero(np.round(t,decimals=3) == np.round(tK[len(tK)-1],decimals=3))
            ende = ende[0][0] + 1
            t      = t[0:ende]
            signal = signal[0:ende]
        elif len(t) < len(tK):
            ende = np.nonzero(np.round(tK,decimals=3) == np.round(tK[len(t)-1],decimals=3))
            ende = ende[0][0] + 1
            tK  = tK[0:ende]
            K   = K[:,0:ende]
    elif factor > 1:
        t      = t[range(0,len(t),factor.astype(np.int32))]
        signal = signal[range(0,len(signal),factor.astype(np.int32))]
        if len(t) > len(tK):
            ende = np.nonzero(np.round(t,decimals=3) == np.round(tK[len(tK)-1],decimals=3))
            ende = ende[0][0] + 1
            t      = t[0:ende]
            signal = signal[0:ende]
        elif len(t) < len(tK):
            ende = np.nonzero(np.round(tK,decimals=3) == np.round(tK[len(t)-1],decimals=3))
            ende = ende[0][0] + 1
            tK  = tK[0:ende]
            K   = K[:,0:ende]
    elif factor < 1:
        tK = tK[range(0,len(tK),(1/factor).astype(np.int32))]
        K  = K[:,range(0,K.shape[1],(1/factor).astype(np.int32))]
        if len(t) > len(tK):
            ende = np.nonzero(np.round(t,decimals=3) == np.round(tK[len(tK)-1],decimals=3))
            ende = ende[0][0] + 1
            t      = t[0:ende]
            signal = signal[0:ende]
        elif len(t) < len(tK):
            ende = np.nonzero(np.round(tK,decimals=3) == np.round(tK[len(t)-1],decimals=3))
            ende = ende[0][0] + 1
            tK  = tK[0:ende]
            K   = K[:,0:ende]

    return t,tK,signal,K

def adjust_nussignal_kernel(t,tK,signal,K):

    lS = len(t)
    lK = len(tK)

    idxK = []
    idxS = []

    if lS < lK:
        for k in range(len(t)):
            idx = np.where(np.round(t[k],3) == np.round(tK,3))
            if len(idx[0]) > 0:
                idxK.append(idx[0][0])
        tK = tK[idxK]
        K = K[:,idxK]
        for k in range(len(tK)):
            idx = np.where(np.round(tK[k],3) == np.round(t,3))
            if len(idx[0]) > 0:
                idxS.append(idx[0][0])
        t = t[idxS]
        signal = signal[idxS]

    elif lS > lK:
        for k in range(len(tK)):
            idx = np.where(np.round(tK[k],3) == np.round(t,3))
            if len(idx[0]) > 0:
                idxS.append(idx[0][0])
        t = t[idxS]
        signal = signal[idxS]
        for k in range(len(tK)):
            idx = np.where(np.round(t[k],3) == np.round(tK,3))
            if len(idx[0]) > 0:
                idxK.append(idx[0][0])
        tK = tK[idxK]
        K = K[:,idxK]

# ------------------------------ #
#     information criterion      #
# ------------------------------ #

def aic(V,Vfit,nrparam):
    N = len(V)
    K = nrparam + 2 # account for variance and intercept of least square fit
    res = np.linalg.norm(V - Vfit)
    value = N*np.log(res**2/N) + 2*K
    return value

def aicc(V,Vfit,nrparam):
    N = len(V)
    K = nrparam + 2 # account for variance and intercept of least square fit
    res = np.linalg.norm(V - Vfit)
    value = N*np.log(res**2/N) + 2*K*N/(N-K-1)
    return value

def bic(V,Vfit,nrparam):
    N = len(V)
    K = nrparam + 2 # account for variance and intercept of least square fit
    res = np.linalg.norm(V - Vfit)
    value = N*np.log(res**2/N) + K*np.log(N)
    return value

def icc(y,fitresult):
    # Get the fitted model
    yfit = fitresult.model
    # Get noise level
    sigma = dl.noiselevel(y[int(np.ceil(len(y)/2)):len(y)])
    # Get non-linear parameters covariance submatrix
    fitpars = fitresult.param + np.finfo(float).eps
    covmat = fitresult.paramUncert.covmat + np.finfo(float).eps
    covmat = covmat/(fitpars[np.newaxis,:]*fitpars[:,np.newaxis])
    # Informational complexity criterion (ICC)
    if not np.all(covmat==0):
        icc = np.sum((y - yfit)**2/sigma**2) + np.sum(np.log(np.diag(covmat))) + np.linalg.slogdet(covmat)[1]
    else:
        icc = np.sum((y - yfit)**2/sigma**2)
    return icc
