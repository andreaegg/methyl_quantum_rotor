import numpy as np
import mqr as mqr
import deerlab as dl

def generatemodel(V3,t,K,bg,dist):

    if bg=='exp':
        if dist=='gauss':
            def fitmodel(Tm,V0,sigma):
                P = mqr.rbd_gauss(V3,V0,sigma)
                B = mqr.bg_exp(t,Tm)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Pfcn = lambda V0,sigma: mqr.rbd_gauss(V3,V0=V0,sigma=sigma)
            Efcn = lambda V0,sigma,scale: K@mqr.rbd_gauss(V3,V0=V0,sigma=sigma)*scale
        elif dist=='2gauss':
            def fitmodel(Tm,V01,sigma1,a1,V02,sigma2,a2):
                P = mqr.rbd_2gauss(V3,V01,sigma1,a1,V02,sigma2,a2)
                B = mqr.bg_exp(t,Tm)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.5, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.5, description='Weight 2')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2: mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,scale: K@mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)*scale
        elif dist=='3gauss':
            def fitmodel(Tm,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
                P = mqr.rbd_3gauss(V3,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3)
                B = mqr.bg_exp(t,Tm)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.3, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.3, description='Weight 2')
            Vmodel.V03.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 3')
            Vmodel.sigma3.set(lb=5, ub=500, par0=100, description='Standard deviation 3')
            Vmodel.a3.set(lb=0, ub=1, par0=0.3, description='Weight 3')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3: mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3,scale: K@mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)*scale
        elif dist=='gengauss':
            def fitmodel(Tm,V0,sigma,beta):
                P = mqr.rbd_gengauss(V3,V0,sigma,beta)
                B = mqr.bg_exp(t,Tm)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.beta.set(lb=0.01, ub=10, par0=2, description='beta parameter')
            Pfcn = lambda V0,sigma,beta: mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)
            Efcn = lambda V0,sigma,beta,scale: K@mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)*scale
        elif dist=='skewgauss':
            def fitmodel(Tm,V0,sigma,alpha):
                P = mqr.rbd_skewgauss(V3,V0,sigma,alpha)
                B = mqr.bg_exp(t,Tm)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.alpha.set(lb=-50, ub=50, par0=0, description='alpha parameter')
            Pfcn = lambda V0,sigma,alpha: mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)
            Efcn = lambda V0,sigma,alpha,scale: K@mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)*scale
        elif dist=='unparam':
            def fitmodel(Tm):
                B = mqr.bg_exp(t,Tm)
                KB = K*B[:,np.newaxis]
                return KB
            Vmodel = dl.Model(fitmodel)
            Vmodel.addlinear('P',vec=len(V3), lb=np.zeros_like(V3))
            Efcn = lambda P: K@P
        else:
            raise KeyError('The chosen distribution model is not available.')
        Vmodel.Tm.set(lb=1, ub=500, par0=5, description='Phase memory time')
        Bfcn = lambda Tm: mqr.bg_exp(t,Tm=Tm)
        
    elif bg=='sumexp':
        if dist=='gauss':
            def fitmodel(Tm1,Tm2,A,V0,sigma):
                P = mqr.rbd_gauss(V3,V0,sigma)
                B = mqr.bg_sumexp(t,Tm1,Tm2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Pfcn = lambda V0,sigma: mqr.rbd_gauss(V3,V0=V0,sigma=sigma)
            Efcn = lambda V0,sigma,scale: K@mqr.rbd_gauss(V3,V0=V0,sigma=sigma)*scale
        elif dist=='2gauss':
            def fitmodel(Tm1,Tm2,A,V01,sigma1,a1,V02,sigma2,a2):
                P = mqr.rbd_2gauss(V3,V01,sigma1,a1,V02,sigma2,a2)
                B = mqr.bg_sumexp(t,Tm1,Tm2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.5, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.5, description='Weight 2')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2: mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,scale: K@mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)*scale
        elif dist=='3gauss':
            def fitmodel(Tm1,Tm2,A,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
                P = mqr.rbd_3gauss(V3,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3)
                B = mqr.bg_sumexp(t,Tm1,Tm2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.3, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.3, description='Weight 2')
            Vmodel.V03.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 3')
            Vmodel.sigma3.set(lb=5, ub=500, par0=100, description='Standard deviation 3')
            Vmodel.a3.set(lb=0, ub=1, par0=0.3, description='Weight 3')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3: mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3,scale: K@mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)*scale
        elif dist=='gengauss':
            def fitmodel(Tm1,Tm2,A,V0,sigma,beta):
                P = mqr.rbd_gengauss(V3,V0,sigma,beta)
                B = mqr.bg_sumexp(t,Tm1,Tm2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.beta.set(lb=0.01, ub=10, par0=2, description='beta parameter')
            Pfcn = lambda V0,sigma,beta: mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)
            Efcn = lambda V0,sigma,beta,scale: K@mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)*scale
        elif dist=='skewgauss':
            def fitmodel(Tm1,Tm2,A,V0,sigma,alpha):
                P = mqr.rbd_skewgauss(V3,V0,sigma,alpha)
                B = mqr.bg_sumexp(t,Tm1,Tm2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.alpha.set(lb=-50, ub=50, par0=0, description='alpha parameter')
            Pfcn = lambda V0,sigma,alpha: mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)
            Efcn = lambda V0,sigma,alpha,scale: K@mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)*scale
        elif dist=='unparam':
            def fitmodel(Tm1,Tm2,A):
                B = mqr.bg_sumexp(t,Tm1,Tm2,A)
                KB = K*B[:,np.newaxis]
                return KB
            Vmodel = dl.Model(fitmodel)
            Vmodel.addlinear('P',vec=len(V3), lb=np.zeros_like(V3))
            Efcn = lambda P: K@P
        else:
            raise KeyError('The chosen distribution model is not available.')
        Vmodel.Tm1.set(lb=1, ub=500, par0=5, description='Phase memory time 1')
        Vmodel.Tm2.set(lb=1, ub=500, par0=5, description='Phase memory time 2')
        Vmodel.A.set(lb=0, ub=1, par0=0.5, description='Weight')
        Bfcn = lambda Tm1,Tm2,A: mqr.bg_sumexp(t,Tm1=Tm1,Tm2=Tm2,A=A)

    elif bg=='prodexp':
        if dist=='gauss':
            def fitmodel(Tm1,Tm2,V0,sigma):
                P = mqr.rbd_gauss(V3,V0,sigma)
                B = mqr.bg_prodexp(t,Tm1,Tm2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Pfcn = lambda V0,sigma: mqr.rbd_gauss(V3,V0=V0,sigma=sigma)
            Efcn = lambda V0,sigma,scale: K@mqr.rbd_gauss(V3,V0=V0,sigma=sigma)*scale
        elif dist=='2gauss':
            def fitmodel(Tm1,Tm2,V01,sigma1,a1,V02,sigma2,a2):
                P = mqr.rbd_2gauss(V3,V01,sigma1,a1,V02,sigma2,a2)
                B = mqr.bg_prodexp(t,Tm1,Tm2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.5, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.5, description='Weight 2')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2: mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,scale: K@mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)*scale
        elif dist=='3gauss':
            def fitmodel(Tm1,Tm2,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
                P = mqr.rbd_3gauss(V3,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3)
                B = mqr.bg_prodexp(t,Tm1,Tm2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.3, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.3, description='Weight 2')
            Vmodel.V03.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 3')
            Vmodel.sigma3.set(lb=5, ub=500, par0=100, description='Standard deviation 3')
            Vmodel.a3.set(lb=0, ub=1, par0=0.3, description='Weight 3')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3: mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3,scale: K@mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)*scale
        elif dist=='gengauss':
            def fitmodel(Tm1,Tm2,V0,sigma,beta):
                P = mqr.rbd_gengauss(V3,V0,sigma,beta)
                B = mqr.bg_prodexp(t,Tm1,Tm2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.beta.set(lb=0.01, ub=10, par0=2, description='beta parameter')
            Pfcn = lambda V0,sigma,beta: mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)
            Efcn = lambda V0,sigma,beta,scale: K@mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)*scale
        elif dist=='skewgauss':
            def fitmodel(Tm1,Tm2,V0,sigma,alpha):
                P = mqr.rbd_skewgauss(V3,V0,sigma,alpha)
                B = mqr.bg_prodexp(t,Tm1,Tm2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.alpha.set(lb=-50, ub=50, par0=0, description='alpha parameter')
            Pfcn = lambda V0,sigma,alpha: mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)
            Efcn = lambda V0,sigma,alpha,scale: K@mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)*scale
        elif dist=='unparam':
            def fitmodel(Tm1,Tm2):
                B = mqr.bg_prodexp(t,Tm1,Tm2)
                KB = K*B[:,np.newaxis]
                return KB
            Vmodel = dl.Model(fitmodel)
            Vmodel.addlinear('P',vec=len(V3), lb=np.zeros_like(V3))
            Efcn = lambda P: K@P
        else:
            raise KeyError('The chosen distribution model is not available.')
        Vmodel.Tm1.set(lb=1, ub=500, par0=5, description='Phase memory time 1')
        Vmodel.Tm2.set(lb=1, ub=500, par0=5, description='Phase memory time 2')
        Bfcn = lambda Tm1,Tm2: mqr.bg_prodexp(t,Tm1=Tm1,Tm2=Tm2)

    elif bg=='strexp':
        if dist=='gauss':
            def fitmodel(Tm,xi,V0,sigma):
                P = mqr.rbd_gauss(V3,V0,sigma)
                B = mqr.bg_strexp(t,Tm,xi)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Pfcn = lambda V0,sigma: mqr.rbd_gauss(V3,V0=V0,sigma=sigma)
            Efcn = lambda V0,sigma,scale: K@mqr.rbd_gauss(V3,V0=V0,sigma=sigma)*scale
        elif dist=='2gauss':
            def fitmodel(Tm,xi,V01,sigma1,a1,V02,sigma2,a2):
                P = mqr.rbd_2gauss(V3,V01,sigma1,a1,V02,sigma2,a2)
                B = mqr.bg_strexp(t,Tm,xi)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.5, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.5, description='Weight 2')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2: mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,scale: K@mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)*scale
        elif dist=='3gauss':
            def fitmodel(Tm,xi,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
                P = mqr.rbd_3gauss(V3,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3)
                B = mqr.bg_strexp(t,Tm,xi)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.3, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.3, description='Weight 2')
            Vmodel.V03.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 3')
            Vmodel.sigma3.set(lb=5, ub=500, par0=100, description='Standard deviation 3')
            Vmodel.a3.set(lb=0, ub=1, par0=0.3, description='Weight 3')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3: mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3,scale: K@mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)*scale
        elif dist=='gengauss':
            def fitmodel(Tm,xi,V0,sigma,beta):
                P = mqr.rbd_gengauss(V3,V0,sigma,beta)
                B = mqr.bg_strexp(t,Tm,xi)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.beta.set(lb=0.01, ub=10, par0=2, description='beta parameter')
            Pfcn = lambda V0,sigma,beta: mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)
            Efcn = lambda V0,sigma,beta,scale: K@mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)*scale
        elif dist=='skewgauss':
            def fitmodel(Tm,xi,V0,sigma,alpha):
                P = mqr.rbd_skewgauss(V3,V0,sigma,alpha)
                B = mqr.bg_strexp(t,Tm,xi)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.alpha.set(lb=-50, ub=50, par0=0, description='alpha parameter')
            Pfcn = lambda V0,sigma,alpha: mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)
            Efcn = lambda V0,sigma,alpha,scale: K@mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)*scale
        elif dist=='unparam':
            def fitmodel(Tm,xi):
                B = mqr.bg_strexp(t,Tm,xi)
                KB = K*B[:,np.newaxis]
                return KB
            Vmodel = dl.Model(fitmodel)
            Vmodel.addlinear('P',vec=len(V3), lb=np.zeros_like(V3))
            Efcn = lambda P: K@P
        else:
            raise KeyError('The chosen distribution model is not available.')
        Vmodel.Tm.set(lb=1, ub=500, par0=5, description='Phase memory time')
        Vmodel.xi.set(lb=0.1, ub=4, par0=2, description='Stretch factor')
        Bfcn = lambda Tm,xi: mqr.bg_strexp(t,Tm=Tm,xi=xi)

    elif bg=='sumstrexp':
        if dist=='gauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,A,V0,sigma):
                P = mqr.rbd_gauss(V3,V0,sigma)
                B = mqr.bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Pfcn = lambda V0,sigma: mqr.rbd_gauss(V3,V0=V0,sigma=sigma)
            Efcn = lambda V0,sigma,scale: K@mqr.rbd_gauss(V3,V0=V0,sigma=sigma)*scale
        elif dist=='2gauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,A,V01,sigma1,a1,V02,sigma2,a2):
                P = mqr.rbd_2gauss(V3,V01,sigma1,a1,V02,sigma2,a2)
                B = mqr.bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.5, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.5, description='Weight 2')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2: mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,scale: K@mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)*scale
        elif dist=='3gauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,A,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
                P = mqr.rbd_3gauss(V3,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3)
                B = mqr.bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.3, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.3, description='Weight 2')
            Vmodel.V03.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 3')
            Vmodel.sigma3.set(lb=5, ub=500, par0=100, description='Standard deviation 3')
            Vmodel.a3.set(lb=0, ub=1, par0=0.3, description='Weight 3')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3: mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3,scale: K@mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)*scale
        elif dist=='gengauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,A,V0,sigma,beta):
                P = mqr.rbd_gengauss(V3,V0,sigma,beta)
                B = mqr.bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.beta.set(lb=0.01, ub=10, par0=2, description='beta parameter')
            Pfcn = lambda V0,sigma,beta: mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)
            Efcn = lambda V0,sigma,beta,scale: K@mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)*scale
        elif dist=='skewgauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,A,V0,sigma,alpha):
                P = mqr.rbd_skewgauss(V3,V0,sigma,alpha)
                B = mqr.bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.alpha.set(lb=-50, ub=50, par0=0, description='alpha parameter')
            Pfcn = lambda V0,sigma,alpha: mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)
            Efcn = lambda V0,sigma,alpha,scale: K@mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)*scale
        elif dist=='unparam':
            def fitmodel(Tm1,xi1,Tm2,xi2,A):
                B = mqr.bg_sumstrexp(t,Tm1,xi1,Tm2,xi2,A)
                KB = K*B[:,np.newaxis]
                return KB
            Vmodel = dl.Model(fitmodel)
            Vmodel.addlinear('P',vec=len(V3), lb=np.zeros_like(V3))
            Efcn = lambda P: K@P
        else:
            raise KeyError('The chosen distribution model is not available.')
        Vmodel.Tm1.set(lb=1, ub=500, par0=5, description='Phase memory time 1')
        Vmodel.xi1.set(lb=0.1, ub=4, par0=2, description='Stretch factor 1')
        Vmodel.Tm2.set(lb=1, ub=500, par0=5, description='Phase memory time 2')
        Vmodel.xi2.set(lb=0.1, ub=4, par0=2, description='Stretch factor 2')
        Vmodel.A.set(lb=0, ub=1, par0=0.5, description='Weight')
        Bfcn = lambda Tm1,xi1,Tm2,xi2,A: mqr.bg_sumstrexp(t,Tm1=Tm1,xi1=xi1,Tm2=Tm2,xi2=xi2,A=A)

    elif bg=='prodstrexp':
        if dist=='gauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,V0,sigma):
                P = mqr.rbd_gauss(V3,V0,sigma)
                B = mqr.bg_prodstrexp(t,Tm1,xi1,Tm2,xi2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Pfcn = lambda V0,sigma: mqr.rbd_gauss(V3,V0=V0,sigma=sigma)
            Efcn = lambda V0,sigma,scale: K@mqr.rbd_gauss(V3,V0=V0,sigma=sigma)*scale
        elif dist=='2gauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,V01,sigma1,a1,V02,sigma2,a2):
                P = mqr.rbd_2gauss(V3,V01,sigma1,a1,V02,sigma2,a2)
                B = mqr.bg_prodstrexp(t,Tm1,xi1,Tm2,xi2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.5, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.5, description='Weight 2')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2: mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,scale: K@mqr.rbd_2gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2)*scale
        elif dist=='3gauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3):
                P = mqr.rbd_3gauss(V3,V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3)
                B = mqr.bg_prodstrexp(t,Tm1,xi1,Tm2,xi2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V01.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 1')
            Vmodel.sigma1.set(lb=5, ub=500, par0=100, description='Standard deviation 1')
            Vmodel.a1.set(lb=0, ub=1, par0=0.3, description='Weight 1')
            Vmodel.V02.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 2')
            Vmodel.sigma2.set(lb=5, ub=500, par0=100, description='Standard deviation 2')
            Vmodel.a2.set(lb=0, ub=1, par0=0.3, description='Weight 2')
            Vmodel.V03.set(lb=1500, ub=2000, par0=1800, description='Distribution mean 3')
            Vmodel.sigma3.set(lb=5, ub=500, par0=100, description='Standard deviation 3')
            Vmodel.a3.set(lb=0, ub=1, par0=0.3, description='Weight 3')
            Pfcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3: mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)
            Efcn = lambda V01,sigma1,a1,V02,sigma2,a2,V03,sigma3,a3,scale: K@mqr.rbd_3gauss(V3,V01=V01,sigma1=sigma1,a1=a1,V02=V02,sigma2=sigma2,a2=a2,V03=V03,sigma3=sigma3,a3=a3)*scale
        elif dist=='gengauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,V0,sigma,beta):
                P = mqr.rbd_gengauss(V3,V0,sigma,beta)
                B = mqr.bg_prodstrexp(t,Tm1,xi1,Tm2,xi2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.beta.set(lb=0.01, ub=10, par0=2, description='beta parameter')
            Pfcn = lambda V0,sigma,beta: mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)
            Efcn = lambda V0,sigma,beta,scale: K@mqr.rbd_gengauss(V3,V0=V0,sigma=sigma,beta=beta)*scale
        elif dist=='skewgauss':
            def fitmodel(Tm1,xi1,Tm2,xi2,V0,sigma,alpha):
                P = mqr.rbd_skewgauss(V3,V0,sigma,alpha)
                B = mqr.bg_prodstrexp(t,Tm1,xi1,Tm2,xi2)
                return (K@P)*B
            Vmodel = dl.Model(fitmodel)
            Vmodel.V0.set(lb=1500, ub=2000, par0=1800, description='Distribution mean')
            Vmodel.sigma.set(lb=5, ub=500, par0=100, description='Standard deviation')
            Vmodel.alpha.set(lb=-50, ub=50, par0=0, description='alpha parameter')
            Pfcn = lambda V0,sigma,alpha: mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)
            Efcn = lambda V0,sigma,alpha,scale: K@mqr.rbd_skewgauss(V3,V0=V0,sigma=sigma,alpha=alpha)*scale
        elif dist=='unparam':
            def fitmodel(Tm1,xi1,Tm2,xi2):
                B = mqr.bg_prodstrexp(t,Tm1,xi1,Tm2,xi2)
                KB = K*B[:,np.newaxis]
                return KB
            Vmodel = dl.Model(fitmodel)
            Vmodel.addlinear('P',vec=len(V3), lb=np.zeros_like(V3))
            Efcn = lambda P: K@P
        else:
            raise KeyError('The chosen distribution model is not available.')
        Vmodel.Tm1.set(lb=1, ub=500, par0=5, description='Phase memory time 1')
        Vmodel.xi1.set(lb=0.1, ub=4, par0=2, description='Stretch factor 1')
        Vmodel.Tm2.set(lb=1, ub=500, par0=5, description='Phase memory time 2')
        Vmodel.xi2.set(lb=0.1, ub=4, par0=2, description='Stretch factor 2')
        Bfcn = lambda Tm1,xi1,Tm2,xi2: mqr.bg_prodstrexp(t,Tm1=Tm1,xi1=xi1,Tm2=Tm2,xi2=xi2)

    else:
        raise KeyError('The chosen background model is not available.')
        
    if dist=='nonparametric':
        return Vmodel,Bfcn,Efcn
    else:
        return Vmodel,Bfcn,Efcn,Pfcn
