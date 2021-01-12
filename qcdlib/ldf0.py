#!/usr/bin/env python
import sys
import os
import numpy as np
import time
from qcdlib.core1 import CORE
#from qcdlib.interpolator import INTERPOLATOR #Either need data files to put in, or create a modified version to handle these
from tools.config import conf

class LDF(CORE):
    """
    upol LDF for electron.
    """

    def __init__(self,lepton='e'):
        self.aux = conf['aux']
        self.set_default_params(lepton)
        self.setup(lepton)
        if lepton=='e': self.ldf=1 #INTERPOLATOR('CJ15lo_0000')
        #if hadron=='pi-': self.pdf=INTERPOLATOR('JAM18PionPDFnlo_0000')

    def set_default_params(self,lepton):

        # free params
        if lepton=='e':
            self._widths1_uv  = 0.
            self._widths1_dv  = 0.
            self._widths1_sea = 0.

            self._widths2_uv  = 0
            self._widths2_dv = 0
            self._widths2_sea = 0

        if lepton=='mu-':
            self._widths1_ubv  = 0.
            self._widths1_dv  = 0.
            self._widths1_sea = 0.

            self._widths2_ubv  = 0
            self._widths2_dv = 0
            self._widths2_sea = 0

        # internal
        self.widths1 = np.ones(11)
        self.widths2 = np.ones(11)

    def setup(self,lepton):
        if lepton=='e':
            for i in range(11):
                if   i == 1: self.widths1[i] = self._widths1_uv
                elif i == 3: self.widths1[i] = self._widths1_dv
                else:        self.widths1[i] = self._widths1_sea
            for i in range(11):
                if   i == 1: self.widths2[i] = self._widths2_uv
                elif i == 3: self.widths2[i] = self._widths2_dv
                else:        self.widths2[i] = self._widths2_sea

        if lepton=='mu-':
            for i in range(11):
                if   i == 2: self.widths1[i] = self._widths1_ubv
                elif i == 3: self.widths1[i] = self._widths1_dv
                else:        self.widths1[i] = self._widths1_sea
            for i in range(11):
                if   i == 2: self.widths2[i] = self._widths2_ubv
                elif i == 3: self.widths2[i] = self._widths2_dv
                else:        self.widths2[i] = self._widths2_sea

    def get_C(self, x, Q2):
        return self.ldf.get_f(x,Q2)

    def get_state(self):
        return self.widths1,self.widths2

    def set_state(self, state):
        self.widths1 = state[0]
        self.widths2 = state[1]

if __name__ == '__main__':

    from qcdlib.aux import AUX

    conf['aux']  = AUX()
    conf['pdfmu-']  = LDF('mu-')
    conf['ldf']  = LDF('e')

    x = 0.15
    Q2 = 2.4
    print conf['ldfmu-'].get_C(x, Q2)
    print conf['ldf'].get_C(x, Q2)
    print conf['aux'].q2qbar(conf['ldf'].get_C(x, Q2))
