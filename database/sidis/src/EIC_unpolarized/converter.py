#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:24:59 2020

@author: avp5627
"""
import os
import pandas as pd
import numpy as np

F=os.listdir('./')
F=[f for f in F if f.endswith('.dat')]

TABLE=[]
for f in F:
  L=open(f).readlines()
  L=[l.strip() for l in L]
  #H=L[1].replace('#PhT','PhT').split()
  H=L[0].split()
  T=L[1:]
  T=[[float(x) for x in l.split()] for l in T if l!='']
  TABLE.extend(T)

PD=pd.DataFrame(TABLE,columns=H)
#print(PD)
PD['sqrts']=np.sqrt(2.*20.*250.)
PD['target']='proton'
PD['obs']='M'
PD['hadron']='pi+'
PD['col']='EIC pseudo'
PD['relative_error_percent']=PD['stat_u']/PD['value']*100.
#print(PD['value'])
PD.to_excel('expdata.xlsx', index=False,header=True)