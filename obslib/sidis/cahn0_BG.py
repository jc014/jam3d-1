#!/usr/bin/env python
import sys
import os
import numpy as np
from tools.tools import load_config
from qcdlib.aux import AUX
from tools.config import conf

eu2, ed2 = 4/9., 1/9.
e2 = []
e2.append(0)    # gluon
e2.append(eu2)  # up
e2.append(eu2)  # ub
e2.append(ed2)  # down
e2.append(ed2)  # db
e2.append(ed2)  # strange
e2.append(ed2)  # sb
e2.append(0)    # charm
e2.append(0)    # cb
e2.append(0)    # bottom
e2.append(0)    # bb
e2 = np.array(e2)

#_get_cahn is used to set parameters and how they vary with charge of particle
def _get_cahn(x, z, Q2, pT, tar, had, F, D, w_tar, w_had):


    if had.endswith('+'):

        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        K = (2*x/(Q2))*(z**2)*(pT**2)*np.abs(w_tar)**2/wq**2
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*F*D*K*gauss)

    elif had.endswith('-'): #change particles from antiparticles to particles

        D=conf['aux'].charge_conj(D)
        w_had=conf['aux'].charge_conj(w_had)
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        K = (2*x/(Q2))*(z**2)*(pT**2)*np.abs(w_tar)**2/wq**2
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*F*D*K*gauss)

    elif had.endswith('0'): #avg of + and - operations

        Dp=D
        Dm=conf['aux'].charge_conj(D)
        w_hadp=w_had
        w_hadm=conf['aux'].charge_conj(w_had)

        wqp = z**2 * np.abs(w_tar) + np.abs(w_hadp)
        wqm = z**2 * np.abs(w_tar) + np.abs(w_hadm)
        gaussp = np.exp(-pT**2 / wqp) / (np.pi * wqp)
        gaussm = np.exp(-pT**2 / wqm) / (np.pi * wqm)
        Kp = (2*x/(Q2))*(z**2)*(pT**2)*np.abs(w_tar)**2/wqp**2
        Km = (2*x/(Q2))*(z**2)*(pT**2)*np.abs(w_tar)**2/wqm**2

        FUUp=np.sum(e2*Kp*F*Dp*gaussp)
        FUUm=np.sum(e2*Km*F*Dm*gaussm)
        return 0.5*(FUUp+FUUm)

def get_cahn(x, z, Q2, pT, tar, had): #uses _get_cahn with the right signs

    # gets collinear parts of positive hadrons and get widths of hadrons
    F = conf['pdf'].get_C(x, Q2)
    w_tar = conf['pdf'].get_widths(Q2)

    if   'pi' in had:
        D = conf['ffpi'].get_C(z, Q2)
        w_had=np.abs(conf['ffpi'].get_widths(Q2))
    elif  'k' in had:
        D = conf['ffk'].get_C(z, Q2)
        w_had=np.abs(conf['ffk'].get_widths(Q2))
    elif 'h' in had:
        D = conf['ffh'].get_C(z,Q2)
        w_had=np.abs(conf['ffh'].get_widths(Q2))
    F[0],D[0]=0,0



    #name function based off target
    if tar == 'p':
        return _get_cahn(x, z, Q2, pT, tar, had, F, D, w_tar, w_had)

    elif tar == 'n':
        F = conf['aux'].p2n(F)
        w_tar = conf['aux'].p2n(w_tar)
        return _get_cahn(x, z, Q2, pT, tar, had, F, D, w_tar, w_had)

    elif tar == 'd':
        return 0.5*(get_cahn(x, z, Q2, pT, 'p', had)+ get_cahn(x, z, Q2, pT, 'n', had))


if __name__ == '__main__':

    from qcdlib.pdf0 import PDF
    from qcdlib.ff0 import FF
    conf['aux']= AUX()
    conf['pdf']=PDF()
    conf['ffpi']=FF('pi')
    conf['ffk']=FF('k')

    x = 0.25
    z = 0.5
    Q2 = 2.4
    mu2 = Q2
    pT = 0.3
    tar = 'p'
    had = 'pi+'

    print get_cahn(x,z,Q2,pT,'p','pi+')
    print get_cahn(x,z,Q2,pT,'p','pi-')
    print get_cahn(x,z,Q2,pT,'p','pi0')

    print get_cahn(x,z,Q2,pT,'n','pi+')
    print get_cahn(x,z,Q2,pT,'n','pi-')
    print get_cahn(x,z,Q2,pT,'n','pi0')

    print get_cahn(x,z,Q2,pT,'d','pi+')
    print get_cahn(x,z,Q2,pT,'d','pi-')
    print get_cahn(x,z,Q2,pT,'d','pi0')
