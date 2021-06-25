import numpy as np
import math as m



def adjust_signal_kernel(t,tK,signal,K):

    factor = (tK[1]-tK[0])/(t[1]-t[0])
    lK = K.shape[1]
    lS = len(signal)

    if lK < lS:
        currt = t[0:lK-1]
        if (t[0] < tK[0]):
            for n in range(lK):
                start = np.nonzero(np.round(currt,decimals=3) == np.round(tK[n],decimals=3))
                if not start:
                    break
            t      = t[start:lS]
            signal = signal[start:lS]
            tK     = tK[n:lK]
            K      = K[:,n:lK]
        elif t[0] > tK[0]:
            for n in range(lK):
                start = np.nonzero(np.round(tK,decimals=3) == np.round(currt[n],decimals=3))
                if not start:
                    break
            tK     = tK[start:lK]
            K      = K[:,start:lK]
            t      = tK[n:lS]
            signal = signal[n:lS]
    
    elif lK > lS:
        currK = tK[0:lS-1]
        if t[0] < tK[0]:
            for n in range(lS):
                start = np.nonzero(np.round(t,decimals=3) == np.round(currK[n],decimals=3))
                if not start:
                    break
            t      = t[start:lS]
            signal = signal[start:lS]
            tK     = tK[n:lK]
            K      = K[:,n:lK]
        elif t[0] > tK[0]:
            for n in range(lS):
                start = np.nonzero(np.round(currK,decimals=3) == np.round(t[n],decimals=3))
                if not start:
                    break
            tK     = tK[start:lK]
            K      = K[:,start:lK]
            t      = t[n:lS-1]
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
