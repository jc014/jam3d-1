#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:00:53 2019

@author: avp5627
"""
import sys
import os
import numpy as np
from tools.tools import load_config
from qcdlib.aux import AUX
from tools.config import conf

eu2, ed2 = 4/9., 1/9.
e2 = []
e2.append(0)    # g
e2.append(eu2)  # u
e2.append(eu2)  # ub
e2.append(ed2)  # d
e2.append(ed2)  # db
e2.append(ed2)  # s
e2.append(ed2)  # sb
e2.append(0)    # c
e2.append(0)    # cb
e2.append(0)    # b
e2.append(0)    # bb
e2 = np.array(e2)

def _get_FUT(xA,xB,Q2,qT,hadronA,hadronB,TransversePolarizationA,TransversePolarizationB,PDFA,PDFB,w_hadronA,w_hadronB):
    """
    We use notations form Phys.Rev. D79 (2009) 034005 http://inspirehep.net/record/796530
    hadronA - moves along +Z in cm frame
    hadronB - moves along -Z in cm frame
    xA - x of a parton in hadronA
    xB - x of a parton in hadronB
    TransversePolarizationA - False True fro ST of hadronA
    TransversePolarizationB - False True fro ST of hadronB
    PDFA - hadronA pdf distributions
    PDFB - hadronB pdf distributions
    w_hadronA - widths of hadronA TMDs
    w_hadronB - widths of hadronB TMDs
    """
 
    MA=conf['aux'].MA
    MB=conf['aux'].MB

    if TransversePolarizationA: # hadronA is transversely polarised
    
            # FTU1 asymmetry equation (95)
            wq = np.abs(w_hadronA) + np.abs(w_hadronB)
            K = -2 * qT * MA / wq
            gauss = np.exp(-qT**2 / wq) / (np.pi * wq)
            return np.sum(e2*K*PDFA*PDFB*gauss)    
    
    elif TransversePolarizationA: # hadronB is transversely polarised
    
            # FUT1 asymmetry equation (98)
            wq = np.abs(w_hadronA) + np.abs(w_hadronB)
            K = 2 * qT * MB / wq
            gauss = np.exp(-qT**2 / wq) / (np.pi * wq)
            return np.sum(e2*K*PDFA*PDFB*gauss)    
    
    else:     
        return 0 # cannot be any asymmetry if none of the particles is polarised
    



def get_FUT(xA,xB,Q2,qT,hadronA,hadronB,TransversePolarizationA,TransversePolarizationB):

    
    if TransversePolarizationA:
        if hadronA == 'p':
            PDFA = (-1.)*conf['sivers'].get_C(xA, Q2) # DY Sivers is opposite to SIDIS
            w_hadronA=conf['sivers'].get_widths(Q2)
        elif hadronA == 'n':
            PDFA = (-1.)*conf['aux'].p2n(conf['sivers'].get_C(xA, Q2)) # DY Sivers is opposite to SIDIS
            w_hadronA=conf['aux'].p2n(conf['sivers'].get_widths(Q2))
        
        if hadronB == 'p':  
            PDFB = conf['pdf'].get_C(xB, Q2)
            w_hadronB = conf['pdf'].get_widths(Q2)
        elif hadronB == 'n':  
            PDFB = conf['aux'].p2n(conf['sivers'].get_C(xB, Q2))
            w_hadronB = conf['aux'].p2n(conf['pdf'].get_widths(Q2))
        elif hadronB == 'pi-':  
            PDFB = conf['aux'].piplus2piminus(conf['pdfpi'].get_C(xB, Q2))
            w_hadronB = conf['aux'].piplus2piminus(conf['pdfpi'].get_widths(Q2))
        elif hadronB == 'pi+':  
            PDFB = conf['pdfpi'].get_C(xB, Q2)
            w_hadronB = conf['pdfpi'].get_widths(Q2) 
            
    elif TransversePolarizationB:
        if hadronB == 'p':
            PDFB = (-1.)*conf['sivers'].get_C(xB, Q2) # DY Sivers is opposite to SIDIS
            w_hadronB=conf['sivers'].get_widths(Q2)
        elif hadronB == 'n':
            PDFA = (-1.)*conf['aux'].p2n(conf['sivers'].get_C(xB, Q2)) # DY Sivers is opposite to SIDIS
            w_hadronA=conf['aux'].p2n(conf['sivers'].get_widths(Q2))
        
        if hadronA == 'p':  
            PDFB = conf['pdf'].get_C(xA, Q2)
            w_hadronB = conf['pdf'].get_widths(Q2)
        elif hadronA == 'n':  
            PDFA = conf['aux'].p2n(conf['sivers'].get_C(xA, Q2))
            w_hadronA = conf['aux'].p2n(conf['pdf'].get_widths(Q2))
        elif hadronA == 'pi-':  
            PDFA = conf['aux'].piplus2piminus(conf['pdfpi'].get_C(xA, Q2))
            w_hadronA = conf['aux'].piplus2piminus(conf['pdfpi'].get_widths(Q2))
        elif hadronA == 'pi+':  
            PDFA = conf['pdfpi'].get_C(xA, Q2)
            w_hadronB = conf['pdfpi'].get_widths(Q2)        
    

    # build structure function
    return _get_FUT(xA,xB,Q2,qT,hadronA,hadronB,TransversePolarizationA,TransversePolarizationB,PDFA,PDFB,w_hadronA,w_hadronB)
    



if __name__ == '__main__':

    from qcdlib.pdf1 import PDF
    conf['aux']= AUX()
    conf['sivers']=PDF()
    conf['pdf']=PDF()
    conf['pdfpi']=PDF('pi')

    xA = 0.25
    xB = 0.5
    Q2 = 16,
    qT = 0.3
    hadronA = 'pi-'
    TransversePolarizationA = False
    hadronB = 'p'
    TransversePolarizationB = True

    print get_FUT(xA,xB,Q2,qT,hadronA,hadronB,TransversePolarizationA,TransversePolarizationB)

