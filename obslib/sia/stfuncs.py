#!/usr/bin/env python
import sys
import os
import numpy as np
from tools.config import conf

# Update from class structure, copy into new file stfuncs0.py

class STFUNCS:

    def __init__(self):
        self.aux = conf['aux']
        eu2, ed2 = 4 / 9., 1 / 9.
        self.e2 = []
        self.e2.append(0)   # g
        self.e2.append(eu2)  # u
        self.e2.append(eu2)  # ub
        self.e2.append(ed2)  # d
        self.e2.append(ed2)  # db
        self.e2.append(ed2)  # s
        self.e2.append(ed2)  # sb
        self.e2.append(eu2)   # c
        self.e2.append(eu2)   # cb
        self.e2.append(ed2)   # b
        self.e2.append(ed2)   # bb
        self.e2 = np.array(self.e2)

        self.Mh = {}
        self.Mh['pi+'] = self.aux.Mpi
        self.Mh['pi-'] = self.aux.Mpi
        self.Mh['k+'] = self.aux.Mk
        self.Mh['k-'] = self.aux.Mk

        self.D = {}
        self.D[1] = {'k1': 'ff', 'k2': 'ff'}
        self.D[2] = {'k1': 'collins', 'k2': 'collins'}

    def get_cc(self, A):
        Acc = np.zeros(A.size)
        Acc[1] = A[2]
        Acc[2] = A[1]
        Acc[3] = A[4]
        Acc[4] = A[3]
        Acc[5] = A[6]
        Acc[6] = A[5]
        Acc[7] = A[8]
        Acc[8] = A[7]
        Acc[9] = A[10]
        Acc[10] = A[9]
        return Acc

    def get_K(self, i, z1, z2, pT, wq, k1, k2, hadron1, hadron2):
        if i == 1:
            if pT != None:
                return np.ones(11)
            else:
                return 1. / (2. * np.pi) * np.ones(11)
        elif i == 2:
            if pT != None:
                return 4 * self.Mh[hadron1]**2 * z1**2 * pT**2 / wq**2
            else:
                return 2. * self.Mh[hadron1]**2 * z1**2 / (np.pi * wq)

    def get_Wq(self, z1, z2, k1, k2, hadron1, hadron2):
        return (z1**2 * self.get_cc(conf[k1].get_widths(1,hadron2))
                + z2**2 * conf[k2].get_widths(1,hadron1)
                ) / z2**2

    def get_gauss(self, z1, z2, pT, k1, k2, wq):
        if pT != None:
            return np.exp(-pT**2 / wq) / (np.pi * wq)
        else:
            return np.ones(11)

    def ZX(self, i, z1, z2, Q2, pT, hadron1, hadron2):
        k1 = self.D[i]['k1']
        k2 = self.D[i]['k2']
        if k1 == None or k2 == None:
            return 0

        D1 = conf[k1].get_C(z1, Q2, hadron1)
        D2 = self.get_cc(conf[k2].get_C(z2, Q2, hadron2))

        Wq = self.get_Wq(z1, z2, k1, k2, hadron1, hadron2)
        gauss = self.get_gauss(z1, z2, pT, k1, k2, Wq)
        K = self.get_K(i, z1, z2, pT, Wq, k1, k2, hadron1, hadron2)
        return np.sum(self.e2 * K * D1 * D2 * gauss)
        # return self.e2[1]*D1[1]*D2[2]*K[1]*gauss[1]\
        #      +self.e2[2]*D1[2]*D2[1]*K[2]*gauss[2]\
        #      +self.e2[3]*D1[3]*D2[4]*K[3]*gauss[3]\
        #      +self.e2[4]*D1[4]*D2[3]*K[4]*gauss[4]\
        #      +self.e2[5]*D1[5]*D2[6]*K[5]*gauss[5]\
        #      +self.e2[6]*D1[6]*D2[5]*K[6]*gauss[6]

if __name__ == '__main__':

    from qcdlib.interpolator import INTERPOLATOR
    from qcdlib.tmdlib import FF
    from qcdlib.tmdlib import COLLINS
    from qcdlib.aux import AUX

    conf['aux']    = AUX()
    conf['cpipff'] = INTERPOLATOR('dsspipNLO_0000')
    conf['cpimff'] = INTERPOLATOR('dsspimNLO_0000')
    conf['cKpff']  = INTERPOLATOR('dssKpNLO_0000')
    conf['cKmff']  = INTERPOLATOR('dssKmNLO_0000')

    conf['lam2'] = 0.4
    conf['Q02']  = 1.0
    conf['ff']      = FF()
    conf['collins'] = COLLINS()

    stfuncs = STFUNCS()
    z1 = 0.5
    z2 = 0.3
    Q2 = 2.5
    pT = 1.0
    hadron1 = 'pi+'
    hadron2 = 'pi-'
    for i in [2]:
        print i, stfuncs.ZX(i, z1, z2, Q2, pT, hadron1, hadron2)
