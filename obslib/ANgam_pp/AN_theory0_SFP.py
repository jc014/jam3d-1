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


#This includes the SFP term to AN in pp -> gam X (see 1410.3448)

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

def get_G(x, Q2): # (f_1T^{\perp(1)}(x)
    return conf['sivers'].get_C(x, Q2)

def get_mandelstam(s, t, u):
# Convenient combinations of the partonic Mandelstam variables
   m['s'] = s
   m['s2'] = s * s
   m['s3'] = s**3.
   m['t'] = t
   m['t2'] = t * t
   m['t3'] = t**3.
   m['u'] = u
   m['u2'] = u * u
   m['u3'] = u**3.
   m['ostu'] = 1. / (s * t * u)
   m['os'] = 1. / s
   m['ot'] = 1. / t
   m['ou'] = 1. / u
   m['st'] = s / t
   m['su'] = s / u
   m['ts'] = t / s
   m['tu'] = t / u
   m['us'] = u / s
   m['ut'] = u / t
   m['st2'] = s**2. / t**2.
   m['su2'] = s**2. / u**2.
   m['ts2'] = t**2. / s**2.
   m['tu2'] = t**2. / u**2.
   m['us2'] = u**2. / s**2.
   m['ut2'] = u**2. / t**2.
   m['os2'] = 1. / s**2.
   m['ot2'] = 1. / t**2.
   m['ou2'] = 1. / u**2.
   m['os3'] = 1. / s**3.
   m['ot3'] = 1. / t**3.
   m['ou3'] = 1. / u**3.
   return m

def get_Hupol(m):
  # Hard parts for the unpolarized cross section
   N_C = 3.0
   Hupol[1] = (2 * (m['s2'] + m['u2']) / (m['t2']*m['u'])) + (2 *m['s'] * (m['u'] - m['s']) / (N_C * m['t']*m['u2']))
   Hupol[2] = (2 * (m['s2'] + m['u2']) / (m['t2']*m['u'])) + ((2 * N_C * m['s'] / (m['u2'])) + (2 * (m['u2'] + m['s']*m['t']) / (N_C * m['s'] * m['t'] * m['u'])))
   Hupol[3] = 2 * (N_C * N_C * m['t'] * m['u'] - (m['s'] * (m['s'] - m['t']))) / ((N_C * N_C - 1) * m['s'] * m['t'] * m['u'])
   Hupol[4] = (2 * (m['s2'] + m['u2']) / (m['t2']*m['u']))
   Hupol[5] = (2 * (m['s2'] + m['u2']) / (m['t2']*m['u']))
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
  Hupol4 = Hupol[4]

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

  G = get_G(x, Q2)
  Gg = G[0]
  Gu = G[1]
  Gub = G[2]
  Gd = G[3]
  Gdb = G[4]
  Gs = G[5]
  Gsb = G[6]
############################################################################
    #uds
  QScs = 0
# a = u
  QScs += ((e2[1] * e2 [2] * Hupol4 *fub * Gu) + (e2[1] * e2 [4] * Hupol4 *fdb * Gu) + (e2[1] * e2 [6] * Hupol4 *fsb * Gu)) + ((e2[1] * e2 [2] * Hupol4 *fu * Gu) + (e2[1] * e2 [4] * Hupol4 *fd * Gu) + (e2[1] * e2 [6] * Hupol4 *fs * Gu)) + (e2[1]**2 * Hupol3 * Gu * fg)
# a = ub
  QScs += ((e2[2] * e2 [1] * Hupol4 *fu * Gub) + (e2[2] * e2 [3] * Hupol4 *fd * Gub) + (e2[2] * e2 [5] * Hupol4 *fs * Gub)) + ((e2[2] * e2 [1] * Hupol4 *fub * Gub) + (e2[2] * e2 [3] * Hupol4 *fdb * Gub) + (e2[2] * e2 [5] * Hupol4 *fsb * Gub)) + (e2[2]**2 * Hupol3 * Gub * fg)
# a = d
  QScs += ((e2[3] * e2 [2] * Hupol4 *fub * Gd) + (e2[3] * e2 [4] * Hupol4 *fdb * Gd) + (e2[3] * e2 [6] * Hupol4 *fsb * Gd)) + ((e2[3] * e2 [2] * Hupol4 *fu * Gd) + (e2[3] * e2 [4] * Hupol4 *fd * Gd) + (e2[3] * e2 [6] * Hupol4 *fs * Gd)) + (e2[3]**2 * Hupol3 * Gd * fg)
# a = db
  QScs += ((e2[4] * e2 [1] * Hupol4 *fu * Gdb) + (e2[4] * e2 [3] * Hupol4 *fd * Gdb) + (e2[4] * e2 [5] * Hupol4 *fs * Gdb)) + ((e2[4] * e2 [1] * Hupol4 *fub * Gdb) + (e2[4] * e2 [3] * Hupol4 *fdb * Gdb) + (e2[4] * e2 [5] * Hupol4 *fsb * Gdb)) + (e2[4]**2 * Hupol3 * Gdb * fg)
# a = s
  QScs += ((e2[5] * e2 [2] * Hupol4 *fub * Gs) + (e2[5] * e2 [4] * Hupol4 *fdb * Gs) + (e2[5] * e2 [6] * Hupol4 *fsb * Gs)) + ((e2[5] * e2 [2] * Hupol4 *fu * Gs) + (e2[5] * e2 [4] * Hupol4 *fd * Gs) + (e2[5] * e2 [6] * Hupol4 *fs * Gs)) + (e2[5]**2 * Hupol3 * Gs * fg)

  QScs += ((e2[6] * e2 [1] * Hupol4 *fu * Gsb) + (e2[6] * e2 [3] * Hupol4 *fd * Gsb) + (e2[6] * e2 [5] * Hupol4 *fs * Gsb)) + ((e2[6] * e2 [1] * Hupol4 *fub * Gsb) + (e2[6] * e2 [3] * Hupol4 *fdb * Gsb) + (e2[6] * e2 [5] * Hupol4 *fsb * Gsb)) + (e2[6]**2 * Hupol3 * Gsb * fg)

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
    den = get_denomint(xF, pT, rs, nx =10)
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
