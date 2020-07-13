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


#AN_theory0.py - program to calculate A_N in pp -> gam X
#This includes the SGP (Qiu-Sterman) term (see 1410.3448)

flavor = ['g','u','ub','d','db','s','sb']
        #  0   1   2    3    4   5    6
target = ['p']

flavdict = {'g': 0, 'u': 1, 'ub': 2,'d': 3, 'db': 4, 's': 5, 'sb': 6}

# Common color factors and fractions
c = {'r3': 1. / 3., 'r4': 0.25, 'r6': 1. / 6., 'r8': 0.125,
         'r9': 1. / 9., 'r18': 1. / 18., 'r24': 1. / 24., 'r27': 1. / 27.}

m = {}
Hupol = {}
HQS = {}
f = {}
ft = {}

eu2, ed2 = 4/9., 1/9.
e2 = []
e2.append(0) #g
e2.append(eu2) #u
e2.append(eu2) #ub
e2.append(ed2) #d
e2.append(ed2) #db
e2.append(ed2) #s
e2.append(ed2) #sb
e2.append(0) #c
e2.append(0) #cb
e2.append(0) #b
e2.append(0) #bb
e2 = np.array(e2)

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

#  @profile
# Calculation of the unpolarized cross section
def get_upolden(x, xF, pT, rs):

  M = conf['aux'].M

  #Defining all of our internal variables at the hadronic and partonic scales#

  if pT > 1.:
    Q = pT
  else:
    Q = 1.

  Q2 = Q * Q
  C_F = 4.0/3.0
  N_C = 3.0
  # Mandelstam variables at the hadron level
  ss = rs**2
  tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) + (xF * ss / 2.0)
  uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) - (xF * ss / 2.0)
  xp = -x * tt / (x * ss + uu)

  # Mandelstam variables at the parton level
  s = x * xp * ss
  t = x * tt
  u = xp * uu

  # Prefactor
  denfac = (1. / (N_C * (x * ss + uu))) * (1. / (x * xp))

  #calling mandelstam variables
  m=get_mandelstam(s, t, u)

  #Calling Hard Factors
  Hupol=get_Hupol(m)
  Hupol1 = Hupol[1]
  Hupol2 = Hupol[2]
  Hupol3 = Hupol[3]

  # Get arrays of the nonperturbative functions
  ft = get_ft(xp, Q2)
  ftg = ft[0]
  ftu = ft[1]
  ftub = ft[2]
  ftd = ft[3]
  ftdb = ft[4]
  fts = ft[5]
  ftsb = ft[6]

  f = get_f(x, Q2)
  fg = f[0]
  fu = f[1]
  fub = f[2]
  fd = f[3]
  fdb = f[4]
  fs = f[5]
  fsb = f[6]
############################################################################

  upol = 0

  upol += ((ftu * fub * 2. * C_F * Hupol1) + (ftub * fg *Hupol3) + (ftg * fu * Hupol2)) * e2[1]

  upol += ((ftub * fu * 2. * C_F * Hupol1) + (ftu * fg *Hupol3) + (ftg * fub * Hupol2)) * e2[2]

  upol += ((ftd * fdb * 2. * C_F * Hupol1) + (ftdb * fg *Hupol3) + (ftg * fd * Hupol2)) * e2[3]

  upol += ((ftdb * fd * 2. * C_F * Hupol1) + (ftd * fg *Hupol3) + (ftg * fdb * Hupol2)) * e2[4]

  upol += ((fts * fsb * 2. * C_F * Hupol1) + (ftsb * fg *Hupol3) + (ftg * fs * Hupol2)) * e2[5]

  upol += ((ftsb * fs * 2. * C_F * Hupol1) + (fts * fg *Hupol3) + (ftg * fsb * Hupol2)) * e2[6]

  return denfac * upol

#  @profile
# Calculation of the fragmentation term in the transversely polarized cross section
def get_polnum(x, xF, pT, rs):

  M = conf['aux'].M

  #Defining all of our internal variables at the hadronic and partonic scales#

  if pT > 1.:
    Q = pT
  else:
    Q = 1.

  Q2 = Q * Q
  C_F = 4.0/3.0
  N_C = 3.0
  # Mandelstam variables at the hadron level
  ss = rs**2
  tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) + (xF * ss / 2.0)
  uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) - (xF * ss / 2.0)
  xp = -x * tt / (x * ss + uu)

  # Mandelstam variables at the parton level
  s = x * xp * ss
  t = x * tt
  u = xp * uu

  # Prefactor
  numfac = (-2.0 * M * pT) * (1. / (x * xp)) / (x * ss + uu)

  #calling mandelstam variables
  m=get_mandelstam(s, t, u)

  #Calling Hard Factors
  Hupol=get_Hupol(m)
  Hupol1 = Hupol[1]
  Hupol2 = Hupol[2]
  Hupol3 = Hupol[3]

  # Get arrays of the nonperturbative functions
  ft = get_ft(xp, Q2)
  ftg = ft[0]
  ftu = ft[1]
  ftub = ft[2]
  ftd = ft[3]
  ftdb = ft[4]
  fts = ft[5]
  ftsb = ft[6]

  f = get_f(x, Q2)
  fg = f[0]
  fu = f[1]
  fub = f[2]
  fd = f[3]
  fdb = f[4]
  fs = f[5]
  fsb = f[6]

  uQS = get_f1Tp(x, Q2)[1]
  ubQS = get_f1Tp(x, Q2)[2]
  dQS = get_f1Tp(x, Q2)[3]
  dbQS = get_f1Tp(x, Q2)[4]
  sQS = get_f1Tp(x, Q2)[5]
  sbQS = get_f1Tp(x, Q2)[6]
############################################################################

  QScs = 0

  QScs += (((-1 / (N_C**2)) * ftub * Hupol1) + ((1 / (2. * C_F)) * ftg *Hupol2)) * (1. / u) * uQS * e2[1]

  QScs += (((-1 / (N_C**2)) * ftu * Hupol1) + ((1 / (2.* C_F)) * ftg *Hupol2)) * (1. / u) * ubQS * e2[2]

  QScs += (((-1 / (N_C**2)) * ftdb * Hupol1) + ((1 / (2. * C_F)) * ftg *Hupol2)) * (1. / u) * dQS * e2[3]

  QScs += (((-1 / (N_C**2)) * ftd * Hupol1) + ((1 / (2.* C_F)) * ftg *Hupol2)) * (1. / u) * dbQS * e2[4]

  QScs += (((-1 / (N_C**2)) * ftsb * Hupol1) + ((1 / (2. * C_F)) * ftg *Hupol2)) * (1. / u) * sQS * e2[5]

  QScs += (((-1 / (N_C**2)) * fts * Hupol1) + ((1 / (2.* C_F)) * ftg *Hupol2)) * (1. / u) * sbQS * e2[6]

  return QScs * numfac


def get_numint(xF, pT, rs, nx = 10):

    C_F = 4.0/3.0
    N_C = 3.0
    # Mandelstam variables at the hadron level
    ss = rs**2
    tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) + (xF * ss / 2.0)
    uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) - (xF * ss / 2.0)

    # Lower limits of the x integration
    xmin = -uu / (ss + tt)

    dnumerdx = np.vectorize(lambda x: get_polnum(x, xF, pT, rs))
    numer = fixed_quad(dnumerdx, xmin, 1., n = nx)[0]
    return numer

def get_denomint(xF, pT, rs, nx = 10):

    C_F = 4.0/3.0
    N_C = 3.0
    # Mandelstam variables at the hadron level
    ss = rs**2
    tt = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) + (xF * ss / 2.0)
    uu = (- rs * np.sqrt( (pT**2) + (xF * xF * ss / 4.0))) - (xF * ss / 2.0)

    # Lower limits of the x integration
    xmin = -uu / (ss + tt)

    ddenomdx = np.vectorize(lambda x: get_upolden(x, xF, pT, rs))
    denom = fixed_quad(ddenomdx, xmin, 1., n = nx)[0]
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
  pT = 2.
  xF = 0.3
  #xF = 2 * pT / rs
  C_F = 4.0/3.0
  N_C = 3.0

  def test():
    num = get_numint(xF, pT, rs, nx = 10)
    den = get_denomint(xF, pT, rs)

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
