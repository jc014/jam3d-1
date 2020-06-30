#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import math
import time
from tools.tools import load_config
from qcdlib.aux import AUX
from tools.config import conf
from scipy.integrate import quad, dblquad, fixed_quad


#AN_theory0.py - program to calculate A_N in pp -> hX
#This includes both the fragmentation term and the QS term (see 1701.09170)

flavor = ['g','u','ub','d','db','s','sb']
        #  0   1   2    3    4   5    6
target = ['p']
hadron = ['pi+','pi-','pi0']

flavdict = {'g': 0, 'u': 1, 'ub': 2,'d': 3, 'db': 4, 's': 5, 'sb': 6}

# Common color factors and fractions
c = {'r3': 1. / 3., 'r4': 0.25, 'r6': 1. / 6., 'r8': 0.125,
         'r9': 1. / 9., 'r18': 1. / 18., 'r24': 1. / 24., 'r27': 1. / 27.}

m = {}
Hupol = {}
HQS = {}
f = {}
ft = {}
d = {}
h = {}
H1p = {}
H = {}

if 'basis' not in conf:
  conf['basis'] = 'default'

def get_f(x, Q2): # Collinear unpolarized PDF
  return conf['pdf'].get_C(x, Q2)

def get_ft(x, Q2): # Collinear unpolarized PDF
  return conf['pdf'].get_C(x, Q2)

def get_f1Tp(x, Q2): # (f_1T^{\perp(1)}(x) - x*df_1T^{\perp(1)}(x)/dx)
    return conf['sivers'].get_C(x, Q2) - x * conf['dsivers'].get_C(x, Q2)

def get_mandelstam(s, t, u):
# Convenient combinations of the partonic Mandelstam variables
   m['st'] = s / t
   m['su'] = s / u
   m['ts'] = t / s
   m['tu'] = t / u
   m['us'] = u / s
   m['ut'] = u / t
   return m

def get_Hupol(m):
  # Hard parts for the unpolarized cross section
   Hupol[1] = m['ut'] + m['tu']
   Hupol[2] = (-m['st']) - m['ts']
   Hupol[3] = (-m['su']) - m['us']
   return Hupol

# def get_HQS(m):
#     fsi = 1. + m['ut']
#     HQS[0] = 0
#     HQS[1] = -m['ou']*(2.+fsi)*(m['st2']+m['ut2'])*c['r18']
#
#     HQS[2] = -m['ou']*(2.-7.*fsi)*(m['su2']+m['tu2'])*c['r18']
#     HQS[3] = -m['ou']*(-10.-fsi)*m['st']*m['su']*c['r27']
#
#     HQS[4] = -m['ou']*(7.+fsi)*(m['st2']+m['ut2'])*c['r18']
#     HQS[5] = -m['ou']*(-1. -7.*fsi)*(m['us2']+m['ts2'])*c['r18']
#
#     HQS[6] = -m['ou']*(-1.-fsi)*m['us']*m['ut']*c['r27']
#
#     HQS[7] = -m['ou']*(7.-2.*fsi)*(m['su2']+m['tu2'])*c['r18']
#     HQS[8] = -m['ou']*(-1. -2.*fsi)*(m['us2']+m['ts2'])*c['r18']
#     HQS[9] = -m['ou']*(-1. -fsi)*m['ts']*m['tu']*c['r27']
#
#     HQS[10] = -m['ou']*c['r6']*c['r9']*(m['tu']+m['ut'])*(1.+18.*m['ts']*m['us'])-m['ou']*fsi*c['r6']*(m['tu']+m['ut'])*(1.-9.*(m['us'])*(m['us']))
#
#     HQS[11] = -m['ou']*c['r4']*c['r4']*(m['su']+m['us'])*(1.-9.*m['ut']*m['ut'])-m['ou']*fsi*c['r8']*c['r18']*(m['su']+m['us'])*(1.+18.*m['st']*m['ut'])
#
#     HQS[12] = -m['ou']*c['r4']*c['r4']*(m['ts']+m['st'])*(1.-9.*m['tu']*m['tu']) + m['ou']*fsi*c['r4']*c['r4']*(m['ts']+m['st'])*(1.-9.*m['su']*m['su'])
#
#     return HQS

#  @profile
# Calculation of the unpolarized cross section
def get_upolden(x, xF, pT, rs):

  M = conf['aux'].M
  Mh = {}
  Mh['pi+'] = conf['aux'].Mpi
  Mh['pi-'] = conf['aux'].Mpi
  Mh['pi0'] = conf['aux'].Mpi
  Mh['k+'] = conf['aux'].Mk
  Mh['k-'] = conf['aux'].Mk
  Mh['jet']=1 #This is a dummy formula so we don't get a runtime error

  if pT > 1.:
    Q = pT
  else:
    Q = 1.

  Q2 = Q * Q
  C_F = 4/3
  N_C = 3
  # Mandelstam variables at the hadron level
  ss = rs**2
  tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) + (xF * ss / 2)
  uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) - (xF * ss / 2)
  x = ((-ss * uu) - (tt * ss) - (tt * uu) - tt) / ((ss**2) + (tt * ss))
  xp = -x * tt / (x * ss + uu)

  # Mandelstam variables at the parton level
  s = x * ss
  t = x * tt
  u = xp * uu

  # Prefactor
  denfac = 1. / (x * ss + uu)

  m=get_mandelstam(s, t, u)
  Hupol=get_Hupol(m)

  Hupol1 = Hupol[1]
  Hupol2 = Hupol[2]
  Hupol3 = Hupol[3]

  # Get arrays of the nonperturbative functions
  #print x,xp
  f = get_f(x, Q2)
  ft = get_ft(xp, Q2)

  fg = f[0]
  fu = f[1]
  fub = f[2]
  fd = f[3]
  fdb = f[4]
  fs = f[5]
  fsb = f[6]

  ftg = ft[0]
  ftu = ft[1]
  ftub = ft[2]
  ftd = ft[3]
  ftdb = ft[4]
  fts = ft[5]
  ftsb = ft[6]

######################
#What about the e**2?
######################
  upol = 0

  upol += ((ftu * fu) * (2 * C_F * Hupol1)) + ((ftu * fg) *Hupol3) + ((ftg * fu) * Hupol2)

  upol += ((ftub * fub) * (2 * C_F * Hupol1)) + ((ftub * fg) *Hupol3) + ((ftg * fub) * Hupol2)

  upol += ((ftd * fd) * (2 * C_F * Hupol1)) + ((ftd * fg) *Hupol3) + ((ftg * fd) * Hupol2)

  upol += ((ftdb * fdb) * (2 * C_F * Hupol1)) + ((ftdb * fg) *Hupol3) + ((ftg * fdb) * Hupol2)

  upol += ((fts * fs) * (2 * C_F * Hupol1)) + ((fts * fg) *Hupol3) + ((ftg * fs) * Hupol2)

  upol += ((ftsb * fsb) * (2 * C_F * Hupol1)) + ((ftsb * fg) *Hupol3) + ((ftg * fsb) * Hupol2)

  return denfac * upol

#  @profile
# Calculation of the fragmentation term in the transversely polarized cross section
def get_polnum(x, xF, pT, rs):

  M = conf['aux'].M
  Mh = {}
  Mh['pi+'] = conf['aux'].Mpi
  Mh['pi-'] = conf['aux'].Mpi
  Mh['pi0'] = conf['aux'].Mpi
  Mh['k+'] = conf['aux'].Mk
  Mh['k-'] = conf['aux'].Mk
  Mh['jet']= 1 #This is a dummy formula so we don't get a runtime error

  Mh = Mh[had]

  if pT > 1.:
    Q = pT
  else:
    Q = 1.

  Q2 = Q * Q
  C_F = 4/3
  N_C = 3
  # Mandelstam variables at the hadron level
  ss = rs**2
  tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) + (xF * ss / 2)
  uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) - (xF * ss / 2)
  x = ((-ss * uu) - (tt * ss) - (tt * uu) - tt) / ((ss**2) + (tt * ss))
  xp = -x * tt / (x * ss + uu)

  # Mandelstam variables at the parton level
  s = x * ss
  t = x * tt
  u = xp * uu

  # Prefactor
  numfac = (2 * M * pT) / (x * ss + uu)

  m=get_mandelstam(s, t, u)
  Hupol=get_Hupol(m)
  
  Hupol1 = Hupol[1]
  Hupol2 = Hupol[2]
  Hupol3 = Hupol[3]

  # Get arrays of the nonperturbative functions
  f = get_f(x, Q2)
  ft = get_ft(xp, Q2)

  ftg = ft[0]
  ftu = ft[1]
  ftub = ft[2]
  ftd = ft[3]
  ftdb = ft[4]
  fts = ft[5]
  ftsb = ft[6]

  fg = f[0]
  fu = f[1]
  fub = f[2]
  fd = f[3]
  fdb = f[4]
  fs = f[5]
  fsb = f[6]

  uQS = (-2./np.pi) * get_f1Tp(x, Q2)[1]
  ubQS = (-2./np.pi) * get_f1Tp(x, Q2)[2]
  dQS = (-2./np.pi) * get_f1Tp(x, Q2)[3]
  dbQS = (-2./np.pi) * get_f1Tp(x, Q2)[4]
  sQS = (-2./np.pi) * get_f1Tp(x, Q2)[5]
  sbQS = (-2./np.pi) * get_f1Tp(x, Q2)[6]

######################
#What about the e**2?
######################
  ffcs = 0

  ffcs += (((-1 / (N_C**2)) * ftu * Hupol1) + ((1 / (2 * C_F)) * ftg *Hupol2)) * (1 / u) * uQS * (np.pi/2.)

  ffcs += (((-1 / (N_C**2)) * ftub * Hupol1) + ((1 / (2 * C_F)) * ftg *Hupol2)) * (1 / u) * ubQS * (np.pi/2.)

  ffcs += (((-1 / (N_C**2)) * ftd * Hupol1) + ((1 / (2 * C_F)) * ftg *Hupol2)) * (1 / u) * dQS * (np.pi/2.)

  ffcs += (((-1 / (N_C**2)) * ftdb * Hupol1) + ((1 / (2 * C_F)) * ftg *Hupol2)) * (1 / u) * dbQS * (np.pi/2.)

  ffcs += (((-1 / (N_C**2)) * fts * Hupol1) + ((1 / (2 * C_F)) * ftg *Hupol2)) * (1 / u) * sQS * (np.pi/2.)

  ffcs += (((-1 / (N_C**2)) * ftsb * Hupol1) + ((1 / (2 * C_F)) * ftg *Hupol2)) * (1 / u) * sbQS * (np.pi/2.)

  if had=='jet': ffcs=0.0

  return ffcs * numfac

def num_integral(x, xp, xF, pT, rs):
    ss = rs**2
    tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) + (xF * ss / 2)
    uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) - (xF * ss / 2)
    x = ((-ss * uu) - (tt * ss) - (tt * uu) - tt) / ((ss**2) + (tt * ss))
    xp = -x * tt / (x * ss + uu)
    return (1 / (x * xp)) * get_polnum(x, xF, pT, rs)
def denom_integral(x, xp, xF, pT, rs):
    ss = rs**2
    tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) + (xF * ss / 2)
    uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) - (xF * ss / 2)
    x = ((-ss * uu) - (tt * ss) - (tt * uu) - tt) / ((ss**2) + (tt * ss))
    xp = -x * tt / (x * ss + uu)
    return (1 / (x * xp)) * get_upolden(x, xF, pT, rs)

def get_numint(xF, pT, rs):

    C_F = 4/3
    N_C = 3
    # Mandelstam variables at the hadron level
    ss = rs**2
    tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) + (xF * ss / 2)
    uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) - (xF * ss / 2)
    x = ((-ss * uu) - (tt * ss) - (tt * uu) - tt) / ((ss**2) + (tt * ss))
    xp = -x * tt / (x * ss + uu)

    # Lower limits of the x integration
    def xmin(uu, tt, ss): return -uu / (ss + tt)

    numer = quad(lambda x: num_integral(x, xp, xF, pT, rs), xmin(uu, tt, ss), 1., limit=100)[0]
    return numer

def get_denomint(xF, pT, rs):

    C_F = 4/3
    N_C = 3
    # Mandelstam variables at the hadron level
    ss = rs**2
    tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) + (xF * ss / 2)
    uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4))) - (xF * ss / 2)
    x = ((-ss * uu) - (tt * ss) - (tt * uu) - tt) / ((ss**2) + (tt * ss))
    xp = -x * tt / (x * ss + uu)

    # Lower limits of the x integration
    def xmin(uu, tt, ss): return -uu / (ss + tt)

    denom = quad(lambda x: denom_integral(x, xp, xF, pT, rs), xmin(uu, tt, ss), 1, limit=100)[0]
    return denom

# def get_vars(rs, n):
#     xF_over_pT = 2 * np.sinh(n) / rs
#     return xF_over_pT
if __name__ == '__main__':


  from qcdlib.ff0 import FF as FF0
  from qcdlib.ff1 import FF as FF1
  from qcdlib.pdf0 import PDF as PDF0
  from qcdlib.pdf1 import PDF as PDF1
  conf['aux']= AUX()
  conf['pdf']=PDF0()
  conf['collinspi']=FF1('pi')
  conf['collinsk']=FF1('k')
  conf['dcollinspi']=FF1('pi','deriv')
  conf['dcollinsk']=FF1('k','deriv')
  conf['Htildepi']=FF1('pi')
  conf['Htildek']=FF1('k')
  conf['transversity']=PDF1()
  conf['sivers']=PDF1()
  conf['dsivers']=PDF1('deriv')
  conf['ffpi']=FF0('pi')
  conf['ffk']=FF0('k')

  rs = 200.
  tar = 'p'
  #had = 'pi+'
  had='jet'
  pT = 2
  xF = 0.2
  N_C = 3

  def test():
    num = get_numint(xF, pT, rs)
    den = get_denomint(xF, pT, rs) / N_C

    AN = num / den
    print AN

  test()

# from timeit import Timer
# t = Timer("test()", "from __main__ import test")
# print 't elapsed ',t.timeit(number=1)

# def test2():
#  den = anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
#  num = anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
#
#  print den,num
#
# from timeit import Timer
# t = Timer("test2()", "from __main__ import test2")
# print 't elapsed ',t.timeit(number=1)

# start = time.time()
# print anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
# print anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
# end = time.time()
# print 'time=',(end-start)

#  start = time.time()
#  test()
#  end = time.time()
#  print 'time=', (end - start)

  # Integration of the numerator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)

  # Integration of the denominator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)
  #den = dblquad(lambda x,z: ANTHEORY().get_dsig(x,z,xF,pT,rs,tar,had),zmin,1.,xmin,lambda x: 1.)

  #AN = num[0]/den[0]
  # print(AN)
