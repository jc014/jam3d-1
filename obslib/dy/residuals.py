#!/usr/bin/env python
import sys
import os
import numpy as np
from tools.residuals import _RESIDUALS
from tools.config import conf
from obslib.dy import upol0 as upol
from obslib.dy import sivers0 as sivers

class RESIDUALS(_RESIDUALS):

    def __init__(self):
        self.reaction = 'dy'
        self.tabs = conf['dy tabs']
        self.setup()

    def _get_theory(self, entry):
        k, i = entry
        xA = self.tabs[k]['xbeam'][i]
        xB = self.tabs[k]['xtarget'][i]

        Q = self.tabs[k]['M'][i]
        Q2 = Q*Q
        qT = self.tabs[k]['qT'][i]

        exp = self.tabs[k]['value'][i]
        hadronB = self.tabs[k]['target'][i]
        TransversePolarizationB=self.tabs[k]['TargetTransversePolarization'][i]
        hadronA = self.tabs[k]['beam'][i]
        TransversePolarizationA=self.tabs[k]['BeamTransversePolarization'][i]
        obs = self.tabs[k]['obs'][i].strip()
        col = self.tabs[k]['col'][i].strip().upper()


        if obs == 'FU1':

            thy = upol.get_FU1(xA,xB,Q2,qT,hadronA,hadronB)


        elif obs == 'AUTsivers':

            # convention factor
            coeff = 1.

            FUT = sivers.get_FUT(xA,xB,Q2,qT,hadronA,hadronB,TransversePolarizationA,TransversePolarizationB)
            FU1 = upol.get_FU1(xA,xB,Q2,qT,hadronA,hadronB)
            thy = coeff * FUT / FU1

        else:
            print 'ERR: exp=%d obs=%s and hadronB=%s not implemented' % (k, obs, hadronB)
            sys.exit()

        return thy

    def gen_report(self, verb=1, level=1):
        """
        verb = 0: Do not print on screen. Only return list of strings
        verv = 1: print on screen the report
        level= 0: only the total chi2s
        level= 1: include point by point
        """

        L = []

        L.append('reaction: %s' % self.reaction)
        # NEEDS SOME WORK HERE

        return L

if __name__ == '__main__':

    from qcdlib import pdf0
    from qcdlib import pdf1
    from qcdlib.aux import AUX
    from reader import READER

    conf['aux']    = AUX()

    conf['pdf']          = pdf0.PDF('p')
    conf['pdfpi-']       = pdf0.PDF('pi-')
    conf['sivers']       = pdf1.PDF()


    conf['datasets']={}
    conf['datasets']['dy']={}

    conf['datasets']['dy']['xlsx']={}

    # COMPASS Sivers
    conf['datasets']['dy']['xlsx'][1000]='dy/expdata/1000.xlsx'  
    conf['datasets']['dy']['xlsx'][1000]='dy/expdata/1001.xlsx'  
    conf['datasets']['dy']['xlsx'][1000]='dy/expdata/1002.xlsx'  
    conf['datasets']['dy']['xlsx'][1000]='dy/expdata/1003.xlsx'  


    conf['datasets']['dy']['norm']={}
    for k in conf['datasets']['dy']['xlsx']: conf['datasets']['dy']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}
    conf['datasets']['dy']['filters']={}

    conf['dy tabs'] = READER().load_data_sets('dy')

    conf['residuals']= RESIDUALS()
    print conf['residuals'].get_residuals()

