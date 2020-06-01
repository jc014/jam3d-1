#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import math
from tools.tools import load_config
from qcdlib.aux import AUX
from tools.config import conf
from scipy.integrate import quad


#AN_theory0.py - program to calculate A_N in ep -> hX (see 1407.5078, 1701.09170)

flavor = ['g','u','ub','d','db','s','sb']
        #  0   1   2    3    4   5    6
target = ['p']
hadron = ['pi+','pi-']

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

def get_H1p(z, Q2, had): # (H_1^{\perp(1)}(z) - z*dH_1^{\perp(1)}(z)/dz)
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


#Code for generic values for AN
#Code to create equation for unpolarized cross-section
def get_unp(z, xF, pT, rs, tar, had):
    S = rs**2
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    Q2 = pT**2
    x = (-1 * U / z) / (S + (T / z))
    s = x*S
    t = (x * T) / z
    u = U / z

    if tar == 'p':
        f = get_f(x, Q2)
    elif tar == 'n':
        f = conf['aux'].p2n(get_f(x, Q2))

    if had.endswith('+'):
        had = had.strip('+')
        d = get_d(z, Q2, had)
    elif had.endswith('-'):
        had = had.strip('-')
        d = conf['aux'].charge_conj(get_d(z, Q2, had))
    elif had.endswith('0'):
        had = had.strip('0')
        plus = get_d(z, Q2, had)
        minus = conf['aux'].charge_conj(get_d(z, Q2, had))
        d = 0.5 * (plus + minus)

    return (1/(z**2)) * (1/x) * (1/(S + (T/z))) * np.sum(e2 * f * d * ((s**2 + u**2) / (t**2)))

#z-integration for unpolarized cross-section
def get_denom(xF, pT, rs, tar, had):
    S = rs*rs
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    zmin = (-1 * (T + U)) / S

    return quad(lambda z: get_unp(z, xF, pT, rs, tar, had), zmin, 1)[0]

#Code to create eqation for polarized cross-section
def get_pol(z, xF, pT, rs, tar, had):
    S = rs**2
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    Q2 = pT**2
    x = (-1 * U / z) / (S + (T / z))
    s = x*S
    t = (x * T) / z
    u = U / z

    if tar == 'p':
        f = get_f(x, Q2)
        f1Tp = get_f1Tp(x, Q2)
        h = get_h(x, Q2)
    elif tar == 'n':
        f = conf['aux'].p2n(get_f(x, Q2))
        f1Tp = conf['aux'].p2n(get_f1Tp(x, Q2))
        h = conf['aux'].p2n(get_h(x, Q2))


    if had.endswith('+'):
        had = had.strip('+')
        d = get_d(z, Q2, had)
        H1p = get_H1p(z, Q2, had)
        H = get_H(z, Q2, had)

    elif had.endswith('-'):
        had = had.strip('-')
        d = conf['aux'].charge_conj(get_d(z, Q2, had))
        H1p = conf['aux'].charge_conj(get_H1p(z, Q2, had))
        H = conf['aux'].charge_conj(get_H(z, Q2, had))

    elif had.endswith('0'):
        had = had.strip('0')
        plusd = get_d(z, Q2, had)
        minusd = conf['aux'].charge_conj(get_d(z, Q2, had))
        d = 0.5 * (plusd + minusd)

        plusH1p = get_H1p(z, Q2, had)
        minusH1p = conf['aux'].charge_conj(get_H1p(z, Q2, had))
        H1p = 0.5 * (plusH1p + minusH1p)

        plusH = get_H(z, Q2, had)
        minusH = conf['aux'].charge_conj(get_H(z, Q2, had))
        H = 0.5 * (plusH + minusH)

    return (1/(z**3)) * (1/x) * ((-1 * 4 * pT) / (S + (T/z))) * np.sum(e2 * ((conf['aux'].M / u) * d * f1Tp * (0.5 * s *(s**2 + u**2) / (t**3)) + (conf['aux'].Mpi / t) * h * (H1p * ((s*u)/(t**2)) + (1/z)*H * ((s/(t**2))*(u-s)))))

#z-integration of polarized cross-section
def get_num(xF, pT, rs, tar, had):
    S = rs*rs
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    zmin = (-1 * (T + U)) / S

    return quad(lambda z: get_pol(z, xF, pT, rs, tar, had), zmin, 1)[0]

#Code for polarized cross-section Qiu-Sterman portion of equation
def get_QS(z, xF, pT, rs, tar, had):
    S = rs**2
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    Q2 = pT**2
    x = (-1 * U / z) / (S + (T / z))
    s = x*S
    t = (x * T) / z
    u = U / z

    if tar == 'p':
        f1Tp = get_f1Tp(x, Q2)
    elif tar == 'n':
        f1Tp = conf['aux'].p2n(get_f1Tp(x, Q2))


    if had.endswith('+'):
        had = had.strip('+')
        d = get_d(z, Q2, had)
    elif had.endswith('-'):
        had = had.strip('-')
        d = conf['aux'].charge_conj(get_d(z, Q2, had))
    elif had.endswith('0'):
        had = had.strip('0')

        plusd = get_d(z, Q2, had)
        minusd = conf['aux'].charge_conj(get_d(z, Q2, had))
        d = 0.5 * (plusd + minusd)
    return (1/(z**3)) * (1/x) * ((-1 * 4 * pT) / (S + (T/z))) * np.sum(e2 * ((conf['aux'].M / u) * d * f1Tp * (0.5 * s *(s**2 + u**2) / (t**3))))

#z-integration of Qiu-Sterman polarized cross-section
def get_numQS(xF, pT, rs, tar, had):
    S = rs*rs
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    zmin = (-1 * (T + U)) / S

    return quad(lambda z: get_QS(z, xF, pT, rs, tar, had), zmin, 1)[0]

#Code for fragmentation of polarized cross-section equation
def get_frag(z, xF, pT, rs, tar, had):
    S = rs**2
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    Q2 = pT**2
    x = (-1 * U / z) / (S + (T / z))
    s = x*S
    t = (x * T) / z
    u = U / z

    if tar == 'p':
        h = get_h(x, Q2)
    elif tar == 'n':
        h = conf['aux'].p2n(get_h(x, Q2))

    if had.endswith('+'):
        had = had.strip('+')
        H1p = get_H1p(z, Q2, had)
        H = get_H(z, Q2, had)
    elif had.endswith('-'):
        had = had.strip('-')
        H1p = conf['aux'].charge_conj(get_H1p(z, Q2, had))
        H = conf['aux'].charge_conj(get_H(z, Q2, had))
    elif had.endswith('0'):
        had = had.strip('0')
        plusH1p = get_H1p(z, Q2, had)
        minusH1p = conf['aux'].charge_conj(get_H1p(z, Q2, had))
        H1p = 0.5 * (plusH1p + minusH1p)

        plusH = get_H(z, Q2, had)
        minusH = conf['aux'].charge_conj(get_H(z, Q2, had))
        H = 0.5 * (plusH + minusH)

    return (1/(z**3)) * (1/x) * ((-1 * 4 * pT) / (S + (T/z))) * np.sum(e2 * ((conf['aux'].Mpi / t) * h * (H1p * ((s*u)/(t**2)) + (1/z)*H * ((s/(t**2))*(u-s)))))

#z-integration of fragmentation portion of polarized cross-section
def get_numfrag(xF, pT, rs, tar, had):
    S = rs*rs
    T = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) + (0.5 * xF * S)
    U = (-1 * rs) * math.sqrt((pT**2) + (xF**2 * 0.25 * S)) - (0.5 * xF * S)
    zmin = (-1 * (T + U)) / S

    return quad(lambda z: get_frag(z, xF, pT, rs, tar, had), zmin, 1)[0]

#Calculation of AN total
def get_AN(xF, pT, rs, tar, had):
    num = get_num(xF, pT, rs, tar, had)
    denom = get_denom(xF, pT, rs, tar, had)

    return num/denom

#Calculation of AN fragmentation
def get_ANfrag(xF, pT, rs, tar, had):
    numfrag = get_numfrag(xF, pT, rs, tar, had)
    denom = get_denom(xF, pT, rs, tar, had)

    return numfrag/denom

#Calculation of AN Qiu-Sterman
def get_ANQS(xF, pT, rs, tar, had):
    numQS = get_numQS(xF, pT, rs, tar, had)
    denom = get_denom(xF, pT, rs, tar, had)

    return numQS/denom


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

  print get_AN(xF, pT, rs, tar, had)
  print get_ANfrag(xF, pT, rs, tar, had)
  print get_ANQS(xF, pT, rs, tar, had)
