from math import log
import numpy as np

class Struct(object):
    pass

def get_advanced_time(t, M_Inf=2.7, R_GW=400.0):
    R = R_GW * (1 + M_Inf/(2*R_GW))**2
    rstar = R + 2*M_Inf*log(R/(2*M_Inf) - 1)
    return t + rstar
def get_retarded_time(t, M_Inf=2.7, R_GW=400.0):
    R = R_GW * (1 + M_Inf/(2*R_GW))**2
    rstar = R + 2*M_Inf*log(R/(2*M_Inf) - 1)
    return t - rstar

def get_chirp_mass(M1, M2):
    M = M1 + M2
    mu = M1*M2/M
    return mu**(3./5.)*M**(2./5.)
def get_masses_from_q_Mchirp(q, Mchirp=1.188):
    assert(np.all(q <= 1))
    M1 = ((1. + q)**(1./5))/(q**(3./5.))*Mchirp
    M2 = q*M1
    return M1, M2

