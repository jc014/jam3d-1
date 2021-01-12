#!/usr/bin/env python
import sys
import os
import numpy as np
import time
from tools.config import conf
#from qcdlib.interpolator import INTERPOLATOR
from qcdlib.core1 import CORE

class FF(CORE):
    """
    upol FF for e-,mu-. Use charge conj 2 to get positive FF
    """

    def __init__(self,lepton):
        self.lepton=lepton
        self.aux = conf['aux']
        self.set_default_params()
        self.setup()
        if   lepton=='e': self.ff = 1. #INTERPOLATOR('dsspipLO_0000')
        elif lepton=='mu':  self.ff = 1. #INTERPOLATOR('dssKpLO_0000')
        #elif lepton=='h':
        #    self.ffpi = INTERPOLATOR('dsspipLO_0000')
        #    self.ffk = INTERPOLATOR('dssKpLO_0000')

    def set_default_params(self):

        # free parameters
        self._widths1_fav   = 0
        self._widths1_ufav  = 0
        self._widths2_fav   = 0
        self._widths2_ufav  = 0

        # internal parameters
        self.widths1 = np.ones(11)
        self.widths2 = np.ones(11)

    def setup(self):
        # 1,  2,  3,  4,  5,  6
        # e, eb,  mu, mub,  tau, taub
        if self.lepton=='mu':
            for i in range(1, 11):
                if   i == 3 or i==4: self.widths1[i] = self._widths1_fav
                else:                  self.widths1[i] = self._widths1_ufav
                if   i == 3 or i==4: self.widths2[i] = self._widths2_fav
                else:                  self.widths2[i] = self._widths2_ufav
        elif self.lepton=='e':
            for i in range(1, 11):
                if   i == 1 or i==2: self.widths1[i] = self._widths1_fav
                else:                  self.widths1[i] = self._widths1_ufav
                if   i == 1 or i==2: self.widths2[i] = self._widths2_fav
                else:                  self.widths2[i] = self._widths2_ufav
        #elif self.lepton=='h':
        #    for i in range(1, 11):
        #        if   i == 1 or i==4 or i==6: self.widths1[i] = self._widths1_fav
        #        else:                  self.widths1[i] = self._widths1_ufav
        #        if   i == 1 or i==4 or i==6: self.widths2[i] = self._widths2_fav
        #        else:                  self.widths2[i] = self._widths2_ufav

    def get_C(self, z, Q2):
        if self.lepton == 'mu' or self.lepton == 'e': return self.ff.get_f(z, Q2)
        #elif self.lepton == 'h': return self.ffpi.get_f(z, Q2) + self.ffk.get_f(z, Q2)

    def get_state(self):
        return self.widths1,self.widths2

    def set_state(self, state):
        self.widths1 = state[0]
        self.widths2 = state[1]

if __name__ == '__main__':

    from qcdlib.aux import AUX
    conf['aux']    = AUX()

    conf['ffmu'] = FF('mu')
    conf['ffe']  = FF('e')
    #conf['ffh']  = FF('h')

    z = 0.15
    Q2 = 2.4
    print conf['ffmu'].get_C(z, Q2)
    print conf['ffe'].get_C(z, Q2)
    #print conf['ffh'].get_C(z, Q2)
