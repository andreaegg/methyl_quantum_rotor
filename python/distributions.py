import numpy as np
import deerlab as dl

def gaussian_V3model(V3,V0,sigma):
    P = dl.dd_gauss(V3,[V0,sigma])
    return P


