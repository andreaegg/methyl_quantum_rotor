import numpy as np
import deerlab as dl
import math as m
import scipy as spy
from scipy.stats import norm

# ------------------------------ #
#      hindered rotor model      #
# ------------------------------ #

def rotorhamiltonian(V3,m,n,B=7.589357490029396):

    if len(m) is not len(n):
        raise(ValueError('Matrix must be symmetric'))

    ham = np.zeros((len(m),len(n)))
    for i in range(len(m)):
        for j in range(len(n)):
            Hij = (n[j]**2)*B*(m[i] == n[j])
            Vij = (V3/4)*(2*(m[i] == n[j]) - ((n[j]-m[i]+3) == 0) - ((n[j]-m[i]+3) == 0))
            ham[i,j] = Hij + Vij

    return ham


def hindered_rotor_model(V3,m,n,B):
    
    E = np.zeros((len(m),len(n)))
    vt = np.zeros((len(V3),1))
    for k in range(len(V3)):
        currpot = V3[k]
        H = rotorhamiltonian(currpot,m,n,B)
        En = np.linalg.eig(H)
        En = np.sort(En)
        E[:,k] = En
        vt[k] = En[1] - En[0]

    return vt,E


def V3pot2vt(V3,nreigenval=49,rotor='CH3'):
    if rotor is 'CH3':
        B = 7.6
    elif rotor is 'CD3':
        B = 3.8

    m = np.arange(-nreigenval,nreigenval+1,1)
    n = np.arange(-nreigenval,nreigenval+1,1)
    
    vt,E = hindered_rotor_model(V3,m,n,B) # vt,E,pot in K 

    return vt,V3,E


def vt2V3pot(vt,nreigenval=49,rotor='CH3',V3start=0):

    if rotor is 'CH3':
        B  = 7.6
    elif rotor is 'CD3':
        B = 3.8
    
    boltzm = 1.3806*10**(-23)
    planck = 6.6261*10**(-34)

    def diff_vtfit(V,vt,m,n,B):
        vtcurr = hindered_rotor_model(V,m,n,B)
        vtcurr = vtcurr*boltzm/planck
        diff = abs(vtcurr-vt)
        return diff

    
    vt = vt*boltzm/planck
    m = np.arange(-nreigenval,nreigenval+1,1)
    n = np.arange(-nreigenval,nreigenval+1,1)
    currV3 = np.zeros((len(vt),1))

    for k in range(len(vt)):
        currvt = vt[k]
        model = lambda x: diff_vtfit(x,currvt,m,n,B)
        V3opt = spy.optimize.fmin(model,V3start)
        currV3[k] = V3opt

    vt,E = hindered_rotor_model(currV3,m,n,B)

    return currV3,vt,E


# ------------------------------ #
# parametric distribution models #
# ------------------------------ #

def normcdf(V,V0,sigma,*args):
    if args:
        alpha = args[0]
    else:
        alpha = 1

    x = (V-V0)/(sigma*np.sqrt(2))
    P = 1/2 + spy.special.erf(alpha*x)/2
    P   = P/np.sum(P)
    return P

def gaussian_V3model(V3,V0,sigma):
    P = dl.dd_gauss(V3,[V0,sigma])
    return P

def rbd_gauss(V,V0,sigma):
    x = (V-V0)/(sigma*np.sqrt(2))
    P = (1/(sigma*np.sqrt(2*m.pi)))*np.exp(-(x**2))
    P   = P/np.trapz(P,V)
    return P

def rbd_2gauss(V,V01,sigma1,a1,V02,sigma2,a2):
    x1 = (V-V01)/(sigma1*np.sqrt(2))
    x2 = (V-V02)/(sigma2*np.sqrt(2))
    P = (a1/(sigma1*np.sqrt(2*m.pi)))*np.exp(-(x1**2)) + (a2/(sigma2*np.sqrt(2*m.pi)))*np.exp(-(x2**2))
    P   = P/np.trapz(P,V)
    return P

def rbd_3gauss(V,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
    x1 = (V-V01)/(sigma1*np.sqrt(2))
    x2 = (V-V02)/(sigma2*np.sqrt(2))
    x3 = (V-V03)/(sigma3*np.sqrt(2))
    P = (a1/(sigma1*np.sqrt(2*m.pi)))*np.exp(-(x1**2)) + (a2/(sigma2*np.sqrt(2*m.pi)))*np.exp(-(x2**2)) + (a3/(sigma3*np.sqrt(2*m.pi)))*np.exp(-(x3**2))
    P   = P/np.trapz(P,V)
    return P

def rbd_gengauss(V,V0,sigma,beta):
    x = abs(V-V0)/sigma
    P = (beta/(2*sigma*spy.special.gamma(1/beta)))*np.exp(-(x**beta))
    P   = P/np.trapz(P,V)
    return P

def rbd_skewgauss(V,V0,sigma,alpha):
    x = (V-V0)/(sigma*np.sqrt(2))
    pdf = (1/(sigma*np.sqrt(2*m.pi)))*np.exp(-(x**2))
    cdf = normcdf(V,V0,sigma,alpha)
    P   = 2*pdf*cdf
    P   = P/np.trapz(P,V)
    return P

def rbd_powergauss(V,V0,sigma,p):
    x = abs(V-V0)/(sigma*np.sqrt(2))
    pdf = (1/(sigma*np.sqrt(2*m.pi)))*np.exp(-(x**2))
    cdf = normcdf(-V,-V0,sigma)
    P   = p*pdf*(cdf**(p-1))
    P   = P/np.trapz(P,V)
    return P

def rbd_lognormal(V,V0,sigma):
    x = (np.log(V)-V0)/(sigma*np.sqrt(2))
    P = (1/(V*sigma*np.sqrt(2*m.pi)))*np.exp(-(x**2))
    P   = P/np.trapz(P,V)
    return P

def rbd_powerlognormal(V,V0,sigma,p):
    x = (np.log(V-V0))/(sigma*np.sqrt(2))
    pdf = (1/(V*sigma*np.sqrt(2*m.pi)))*np.exp(-(x**2))
    cdf = normcdf(-np.log(V),-np.log(V0),sigma)
    P   = p*pdf*(cdf**(p-1))
    P   = P/np.trapz(P,V)
    return P


# ------------------------------ #
#        background models       #
# ------------------------------ #

def bg_exp(t,Tm):
    bg = np.exp(-t/Tm)
    return bg

def bg_strexp(t,Tm,xi):
    bg = np.exp(-((t/Tm)**xi))
    return bg

def bg_sumexp(t,Tm1,Tm2,A):
    bg = A*np.exp(-t/Tm1) + (1-A)*np.exp(-t/Tm2)
    return bg

def bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A):
    bg = A*np.exp(-((t/Tm1)**xi1)) + (1-A)*np.exp(-((t/Tm2)**xi2))
    return bg    

def bg_prodexp(t,Tm1,Tm2):
    bg = np.exp(-t/Tm1)*np.exp(-t/Tm2)
    return bg

def bg_prodstrexp(t,Tm1,xi1,Tm2,xi2):
    bg = np.exp(-((t/Tm1)**xi1))*np.exp(-((t/Tm2)**xi2))
    return bg  

def bg_2ndpolynomial(t,a0,a1,a2):
    bg = a0 + a1*(t**(-1)) + a2*(t**(-2))
    return bg

def bg_3rdpolynomial(t,a0,a1,a2,a3):
    bg = a0 + a1*(t**(-1)) + a2*(t**(-2)) + a3*(t**(-3))
    return bg

def bg_4thpolynomial(t,a0,a1,a2,a3,a4):
    bg = a0 + a1*(t**(-1)) + a2*(t**(-2)) + a3*(t**(-3)) + a4*(t**(-4))
    return bg

def bg_5thpolynomial(t,a0,a1,a2,a3,a4,a5):
    bg = a0 + a1*(t**(-1)) + a2*(t**(-2)) + a3*(t**(-3)) + a4*(t**(-4)) + a5*(t**(-5))
    return bg

def bg_6thpolynomial(t,a0,a1,a2,a3,a4,a5,a6):
    bg = a0 + a1*(t**(-1)) + a2*(t**(-2)) + a3*(t**(-3)) + a4*(t**(-4)) + a5*(t**(-5)) + a6*(t**(-6))
    return bg