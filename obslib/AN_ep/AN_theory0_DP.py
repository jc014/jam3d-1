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


#AN_theory0_DP.py - cross-check of program to calculate A_N in ep -> hX (see 1407.5078, 1701.09170)

flavor = ['g','u','ub','d','db','s','sb']
        #  0   1   2    3    4   5    6
target = ['p']
hadron = ['pi+','pi-','pi0']

flavdict = {'g': 0, 'u': 1, 'ub': 2,'d': 3, 'db': 4, 's': 5, 'sb': 6}

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

f = {}
d = {}
h = {}
H1p = {}
H = {}
f1Tp = {}

if 'basis' not in conf:
  conf['basis'] = 'default'

def get_f(x, Q2): # Collinear unpolarized PDF
  return conf['pdf'].get_C(x, Q2)

def get_d(z, Q2, had): # Collinear unpolarized FF
  if 'pi' in had:
      return conf['ffpi'].get_C(z, Q2)
  elif 'k' in had:
      return conf['ffk'].get_C(z, Q2)

def get_h(x, Q2): # Collinear transversity
  return conf['transversity'].get_C(x, Q2)

# (H_1^{\perp(1)}(z) - z*dH_1^{\perp(1)}(z)/dz)
def get_H1p(z, Q2, had):
  if 'pi' in had:
      return conf['collinspi'].get_C(z, Q2) - z * conf['dcollinspi'].get_C(z, Q2)
  elif 'k' in had:
      return conf['collinsk'].get_C(z, Q2) - z * conf['dcollinsk'].get_C(z, Q2)

def get_H(z, Q2, had): # -2*z*H_1^{\perp(1)}(z)+\tilde{H}(z)
  if 'pi' in had:
      return -2. * z * conf['collinspi'].get_C(z,Q2) + conf['Htildepi'].get_C(z, Q2)
  elif 'k' in had:
      return -2. * z * conf['collinsk'].get_C(z,Q2) + conf['Htildek'].get_C(z, Q2)

def get_f1Tp(x, Q2): # (f_1T^{\perp(1)}(x) - x*df_1T^{\perp(1)}(x)/dx)
    return conf['sivers'].get_C(x, Q2) - x * conf['dsivers'].get_C(x, Q2)

def get_Hupol(s, t, u):
  # Hard parts for the unpolarized cross section
   Hupol = (s**2 + u**2)/t**2
   return Hupol

def get_HTffa(s, t, u):
  # Hard parts for the transversely polarized fragmentation term
   HTffa = s*u/t**2
   return HTffa

def get_HTffb(s, t, u):
  # Hard parts for the transversely polarized fragmentation term
   HTffb = s*(u-s)/t**2
   return HTffb

def get_Hxxpz(z, Q2, had, s, t, u):
  HTffa = get_HTffa(s, t, u)
  HTffb = get_HTffb(s, t, u)

  H1p = get_H1p(z, Q2, 'pi+')
  H = get_H(z, Q2, 'pi+')

  if had.endswith('-'):
      H1p = conf['aux'].charge_conj(H1p)
      H = conf['aux'].charge_conj(H)

  elif had.endswith('0'):
      H1pp=H1p
      H1pm=conf['aux'].charge_conj(H1p)
      Hp=H
      Hm=conf['aux'].charge_conj(H)

      H1p=0.5*(H1pp+H1pm)
      #print z,2.0*z**2*0.135*H1pp,2.0*z**2*0.135*H1pm,2.0*z**2*0.135*H1p
      #sys.exit()
      H=0.5*(Hp+Hm)

  #print had,z,H1p[1],H1p[3]

  Hxxpz = HTffa * H1p + HTffb * H / z
  return Hxxpz

def get_HQS(s,t,u):
    HQS = s * (s**2 + u**2)/(2.*t**3)
    return HQS

#  @profile
# Calculation of the unpolarized cross section
def get_dsig(z, xF, pT, rs, tar, had):

  M = conf['aux'].M
  Mh = {}
  Mh['pi+'] = conf['aux'].Mpi
  Mh['pi-'] = conf['aux'].Mpi
  Mh['pi0'] = conf['aux'].Mpi
  Mh['k+'] = conf['aux'].Mk
  Mh['k-'] = conf['aux'].Mk


  if pT > 1.:
    Q = pT
  else:
    Q = 1.

  Q2 = Q * Q

  xT = 2. * pT / rs
  xT2 = xT * xT
  xF2 = xF * xF

  # Mandelstam variables at the hadron level
  ss = rs * rs
  tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
  uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

  oz = 1. / z
  x = -(uu/z) / (ss + tt/z)

  # Mandelstam variables at the parton level
  s = x * ss
  t = x * tt * oz
  u = uu * oz

  # Prefactor
  denfac = 1. / (z**2 * x * (ss + tt/z))

  Hupol=get_Hupol(s,t,u)

  # Get arrays of the nonperturbative functions
  f = get_f(x, Q2)
  d = get_d(z, Q2, 'pi+')

  if tar == 'n':
      f = conf['aux'].p2n(f)

  if had.endswith('-'):
      d = conf['aux'].charge_conj(d)

  elif had.endswith('0'):
      dp=d
      dm=conf['aux'].charge_conj(d)
      d=0.5*(dp+dm)
      #print z,z*dp,z*dm,z*d

  upol = np.sum(e2 * f * d * Hupol)

  return denfac * upol

#  @profile
# Calculation of the fragmentation term in the transversely polarized cross section
def get_dsigST(z, xF, pT, rs, tar, had):

  M = conf['aux'].M
  Mh = {}
  Mh['pi+'] = conf['aux'].Mpi
  Mh['pi-'] = conf['aux'].Mpi
  Mh['pi0'] = conf['aux'].Mpi
  Mh['k+'] = conf['aux'].Mk
  Mh['k-'] = conf['aux'].Mk

  Mh = Mh[had]

  if pT > 1.:
    Q = pT
  else:
    Q = 1.

  Q2 = Q * Q

  xT = 2. * pT / rs
  xT2 = xT * xT
  xF2 = xF * xF

  # Mandelstam variables at the hadron level
  ss = rs * rs
  tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
  uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

  oz = 1. / z
  x = -(uu/z) / (ss + tt/z)

  # Mandelstam variables at the parton level
  s = x * ss
  t = x * tt * oz
  u = uu * oz

  # Prefactor
  numfac = -4 * pT / (z**3 * x * (ss + tt/z))

  # Get arrays of the nonperturbative functions
  f = get_f(x, Q2)
  h = get_h(x, Q2)
  f1Tp = get_f1Tp(x, Q2)
  Hxxpz = get_Hxxpz(z, Q2, had, s, t, u)
  HQS = get_HQS(s,t,u)
  d = get_d(z, Q2, 'pi+')

  if tar == 'n':
      f = conf['aux'].p2n(f)
      h = conf['aux'].p2n(h)
      f1Tp = conf['aux'].p2n(f1Tp)

  if had.endswith('-'):
      d = conf['aux'].charge_conj(d)

  elif had.endswith('0'):
      dp=d
      dm=conf['aux'].charge_conj(d)
      d=0.5*(dp+dm)

  ffcs = np.sum(e2 * (Mh/t) * h * Hxxpz)


  #Qiu-Sterman term

  QScs = np.sum(e2 * (M/u) * d * f1Tp * HQS)

  return numfac * (ffcs + QScs)
  #return ffcs
  #return  QScs

def get_sig(xF, pT, rs, tar, had):
    xT = 2. * pT / rs
    xF2 = xF * xF
    xT2 = xT * xT

    # Mandelstam variables at the hadron level
    ss = rs * rs
    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    # Lower limits of the z and x integrations
    zmin = np.sqrt(xF2 + xT2)

    sig = quad(lambda z: get_dsig(z, xF, pT, rs, tar, had), zmin, 1)[0]

#    if mode == 'gauss':
#        dsigdzdx = np.vectorize(
#            lambda x, z: get_dsig(x, z, xF, pT, rs, tar, had))
#        dsigdz = np.vectorize(lambda z: fixed_quad(
#            lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
#        sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
#    elif mode == 'quad':
#        sig = dblquad(lambda x, z: get_dsig(
#            x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]
    return sig

def get_sigST(xF, pT, rs, tar, had):
    xT = 2. * pT / rs
    xF2 = xF * xF
    xT2 = xT * xT

    # Mandelstam variables at the hadron level
    ss = rs * rs
    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    # Lower limits of the z and x integrations
    zmin = np.sqrt(xF2 + xT2)

    sig = quad(lambda z: get_dsigST(z, xF, pT, rs, tar, had), zmin, 1)[0]

#    if mode == 'gauss':
#      dsigdzdx = np.vectorize(
#        lambda x, z: get_dsigST(x, z, xF, pT, rs, tar, had))
#      dsigdz = np.vectorize(lambda z: fixed_quad(
#        lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
#      sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
#    elif mode == 'quad':
#      sig = dblquad(lambda x, z: get_dsigST(
#        x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]

    #print had,sig
    return sig


if __name__ == '__main__':

    from qcdlib.ff0 import FF as FF0
    from qcdlib.ff1 import FF as FF1
    from qcdlib.pdf0 import PDF as PDF0
    from qcdlib.pdf1 import PDF as PDF1
    conf['aux']= AUX()
    conf['pdf']=PDF0()
    conf['collinspi']=FF1('pi')
    conf['dcollinspi']=FF1('pi','deriv')
    conf['dcollinsk']=FF1('k','deriv')
    conf['Htildepi']=FF1('pi')
    conf['transversity']=PDF1()
    conf['sivers']=PDF1()
    conf['dsivers']=PDF1('deriv')
    conf['ffpi']=FF0('pi')

    rs = 7.25
    tar = 'p'
    had = 'pi+'
    pT = 1.0
    xF = 0.3
    Q2 = 2.0

    num = get_sigST(xF, pT, rs, tar, had)
    den = get_sig(xF, pT, rs, tar, had)
    AN = num/den

    print AN


#  from qcdlib.ff0 import FF as FF0
#  from qcdlib.ff1 import FF as FF1
#  from qcdlib.pdf0 import PDF as PDF0
#  from qcdlib.pdf1 import PDF as PDF1
#  conf['aux']= AUX()
#  conf['pdf']=PDF0()
#  conf['collinspi']=FF1('pi')
#  conf['collinsk']=FF1('k')
#  conf['Htildepi']=FF1('pi')
#  conf['Htildek']=FF1('k')
#  conf['transversity']=PDF1()
#  conf['sivers']=PDF1()
#  conf['ffpi']=FF0('pi')
#  conf['ffk']=FF0('k')
#
#  rs = 200.
#  tar = 'p'
#  had = 'pi+'
#  pT = 3.0
#  xF = 0.2375
#  xT = 2. * pT / rs
#  xF2 = xF * xF
#  xT2 = xT * xT
#
#  # Mandelstam variables at the hadron level
#  ss = rs * rs
#  tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
#  uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)
#
#  # Lower limits of the z and x integrations
#  zmin = np.sqrt(xF2 + xT2)
#
#  def xmin(z): return -uu / (z * ss + tt)
#
#  def test():
#    den = get_sig(xF, pT, rs, tar, had, mode='gauss', nx=100, nz=100)
#    num = get_sigST(xF, pT, rs, tar, had,
#               mode='gauss', nx=100, nz=100)
#
#    AN = num / den
#    print AN
#
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
