#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#change all Mandelstam definitons, etc.
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

#ps_RC: calculate parton scattering unpolarized cross-section

flavor = ['g','u','ub','d','db','s','sb']
    #  0  1  2  3  4  5  6
target = ['p']
hadron = ['pi+','pi-','pi0']

flavdict = {'g': 0, 'u': 1, 'ub': 2,'d': 3, 'db': 4, 's': 5, 'sb': 6}

# Common color factors and fractions
c = {'r3': 1. / 3., 'r4': 0.25, 'r6': 1. / 6., 'r8': 0.125,
         'r9': 1. / 9., 'r18': 1. / 18., 'r24': 1. / 24., 'r27': 1. / 27.}

#charges
e=[0,2/3,-2/3,-1/3,1/3,-1/3,1/3]

m = {}
rf = {}
Hupol = {}
f = {}
ft = {}
d = {}
h = {}
H1p = {}
H = {}
F1p = {}
#HTffa = np.zeros(13)
#HTffb = np.zeros(13)
#Hxxpz = np.zeros((13, 7))

if 'basis' not in conf:
    conf['basis'] = 'default'
#need to add options for LDF, LFF
def get_f(x, Q2): # Collinear unpolarized PDF
  return conf['pdf'].get_C(x, Q2)

def get_ft(x, Q2): # Collinear unpolarized PDF
  return conf['pdf'].get_C(x, Q2)

def get_d(z, Q2, had): # Collinear unpolarized FF
  if 'pi' in had:
      return conf['ffpi'].get_C(z, Q2)
  elif 'k' in had:
      return conf['ffk'].get_C(z, Q2)

def get_h(x, Q2): # Collinear transversity
  return conf['transversity'].get_C(x, Q2)

def get_fl(xi, Q2):
    #return 1.
    return conf['ldf'].get_C(xi, Q2)

def get_dl(zt, Q2):
    #return 1.
    return conf['lff'].get_C(zt, Q2)

## (H_1^{\perp(1)}(z) - z*dH_1^{\perp(1)}(z)/dz)
#def get_H1p(z, Q2, had):
#  if 'pi' in had:
#      return conf['collinspi'].get_C(z, Q2) - z * conf['collinspi'].get_dC(z, Q2)
#  elif 'k' in had:
#      return conf['collinsk'].get_C(z, Q2) - z * conf['collinsk'].get_dC(z, Q2)
#
#def get_H(z, Q2, had): # -2*z*H_1^{\perp(1)}(z)+\tilde{H}(z)
#  if 'pi' in had:
#      return -2. * z * conf['collinspi'].get_C(z,Q2) + conf['Htildepi'].get_C(z, Q2)
#  elif 'k' in had:
#      return -2. * z * conf['collinsk'].get_C(z,Q2) + conf['Htildek'].get_C(z, Q2)
#
#
## (f_1T^{\perp(1)}(x) - x*df_1T^{\perp(1)}(x)/dx) Check this for correct Sivers implementation
#def get_f1Tp( x, Q2):
#    return conf['sivers'].get_C(x, Q2) - x * conf['sivers'].get_dC(x, Q2)

def get_mandelstam(s, t, u, sp, tp, up):
# Convenient combinations of the partonic Mandelstam variables
   m['s2'] = s * s
   m['s3'] = s**3.
   m['t2'] = t * t
   m['t3'] = t**3.
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
   m['sp2'] = sp * sp
   m['sp3'] = sp**3.
   m['tp2'] = tp * tp
   m['tp3'] = tp**3.
   m['up2'] = up * up
   m['up3'] = up**3.
   m['osptpup'] = 1. / (sp * tp * up)
   m['osp'] = 1. / sp
   m['otp'] = 1. / tp
   m['oup'] = 1. / up
   m['sptp'] = sp / tp
   m['spup'] = sp / up
   m['tpsp'] = tp / sp
   m['tpup'] = tp / up
   m['upsp'] = up / sp
   m['uptp'] = up / tp
   m['sptp2'] = sp**2. / tp**2.
   m['spup2'] = sp**2. / up**2.
   m['tpsp2'] = tp**2. / sp**2.
   m['tpup2'] = tp**2. / up**2.
   m['upsp2'] = up**2. / sp**2.
   m['uptp2'] = up**2. / tp**2.
   m['osp2'] = 1. / sp**2.
   m['otp2'] = 1. / tp**2.
   m['oup2'] = 1. / up**2.
   m['osp3'] = 1. / sp**3.
   m['otp3'] = 1. / tp**3.
   m['oup3'] = 1. / up**3.
   return m

def get_rapidityf(rs,ppT,lpT,yL,yP,thLP):
    rf['ttp']=-rs*ppT*math.exp(yP)
    rf['uup']=-rs*lpT*math.exp(yL)
    rf['tt']=rs*lpT*math.exp(-yL)
    rf['uup']=rs*ppT*math.exp(-yP)
    rf['ss']=rs*rs
    rf['ssp']=ppT*lpT*(math.exp(yL-yP)+math.exp(yP-yL)-2*math.cos(thLP))
    return rf

def get_Hupol(m,pT,zt,xi,x,z,s,t,u,sp,tp,up):
  # Hard parts for the unpolarized cross section
   if pT > 1.:
     Q = pT
   else:
     Q = 1.

   Q2 = Q * Q


   Hupol[1] = 2*Q2*Q2*z*z
   Hupol[2] = 2*zt*Q2*up*z
   Hupol[3] = zt*zt*(m['s2']*x*x*z*z+2*m['up2'])
   Hupol[4] = 2*xi*u*x*z*(Q2*z-zt*up)
   Hupol[5] = m['u2']*x*x*z*z
   return Hupol


# Calculation of the unpolarized cross section
def get_dsig(ppT,lpT,yL,yP,thLP,x, z, xi, zt, pT, rs, tar, had):

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

  # Mandelstam variables at the hadron level
  ss = rs * rs
  tt = -Q2
  #uu =

  oz = 1. / z
  ozt= 1. / zt

  get_rapidityf(rs,ppT,lpT,yL,yP,thLP)

  # Mandelstam variables at the parton level
  s = x * xi * rf['ss']
  t = xi*ozt*rf['tt']
  u = x * rf['uu'] * ozt
  sp = oz*ozt*rf['ssp']
  tp = x*oz*rf['ttp']
  up = xi*oz*rf['uup']

  # Prefactor
  denfacint = 1. / (zt*zt*xi*z*z*x)
  pifac = 2.*math.pi * (1. / (2.*(2.*math.pi)**3.))**2.
  denfac = 1. / (-xi*xi*rf['tt']*rf['ss']*rf['ss']*rf['uu']*x*x*z*z)
  #numfac=-2*el*el*g*g    -fill in constants?
  prefac=denfacint*pifac*denfac

  m=get_mandelstam(s, t, u,sp,tp,up)
  get_Hupol(m,pT,zt,xi,x,z,s,t,u,sp,tp,up)

  Hupol1 = Hupol[1]
  Hupol2 = Hupol[2]
  Hupol3 = Hupol[3]
  Hupol4 = Hupol[4]
  Hupol5 = Hupol[5]

  # Get arrays of the nonperturbative functions
  f = get_f(x, Q2)
  #ft = get_ft(xp, Q2)
  d = get_d(z, Q2, 'pi+')

  fl = get_fl(xi,Q2)

  dl = get_dl(zt,Q2)

  if had.endswith('-'):
      d = conf['aux'].charge_conj(d)

  elif had.endswith('0'):
      dp=d
      dm=conf['aux'].charge_conj(d)
      d=0.5*(dp+dm)
      #print z,z*dp,z*dm,z*d

  fg = f[0]
  fu = f[1]
  fub = f[2]
  fd = f[3]
  fdb = f[4]
  fs = f[5]
  fsb = f[6]


  dg = d[0]
  du = d[1]
  dub = d[2]
  dd = d[3]
  ddb = d[4]
  ds = d[5]
  dsb = d[6]
#e[i] is charge of pdf function
  Hadprod = 0
    for i in range (6):
        for j in range ( i+1,6):
            Hadprod = Hadprod +e[i]*el*e[i]*el*f[i]*d[j]

  Hupol = xi*xi*(Hupol1-Hupol2+Hupol3)+Hupol4+Hupol5

#add ff and pdf terms, and ldf, lff terms (all possible multiplications from sums, with charges)

  upol=fl*dl*Hupol*Hadprod

  return prefac * upol


#adjust this to  4 integrations
def get_sig(lpT,ppT,yL,yP,thLP, pT, rs, tar, had, mode='gauss', nx=10, nz=10,nxi=10,nzt=10):
    #xT = 2. * pT / rs
    #xF2 = xF * xF
    #xT2 = xT * xT

    # Mandelstam variables at the hadron level
    #ss = rs * rs
    #tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    #uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    get_rapidityf(rs,ppT,lpT,yL,yP,thLP)

    # Lower limits of the  integrations
    ztmin = (-rf['tt']-rf['ssp']-rf['uup'])/(rf['ttp']+rf['ss']+rf['uu'])

    def ximin(zt): return (rf['ttp']+(1./zt)*rf['ssp']+(1./zt)*rf['uup'])/((1./zt)*rf['t']-rf'[ss']-rf['uu'])

    def zmin(xi,zt): return (rf['ttp']+(1/zt)*rf['ssp']+xi*rf['uu'])/((xi/zt)*rf['t']-(1/zt)*rf['uup']-xi*rf['ss'])
#x fixed from delta function??
    def xmin(z,xi,zt): return ((xi/zt)*rf['t']-(1/(z*zt))*rf['ssp']+(xi/z)*rf['uu'])/((1/zt)*rf['uup']+xi*rf['ss']+(1/z)*rf['ttp'])

    if mode == 'gauss':
        dsigdztdxidzdx = np.vectorize(
            lambda x, z, xi, zt: get_dsig(ppT,lpT,yL,yP,thLP,x, z, xi,zt, pT, rs, tar, had))
        dsigdztdxidz = np.vectorize(lambda z, xi, zt: fixed_quad(
            lambda x: dsigdztdxidzdx(x, z, xi, zt), xmin(z,xi,zt), 1, n=nx)[0])
        dsigdztdxi = np.vectorize(lambda xi, zt: fixed_quad(
            lambda z: dsigdztdxidzdx(z, xi, zt), zmin(xi,zt), 1, n=nx)[0])
        dsigdzt = np.vectorize(lambda xi: fixed_quad(
            lambda xi: dsigdztdxidz(xi, zt), ximin(zt), 1, n=nxi)[0])
        sig = fixed_quad(dsigdzt, ztmin, 1, n=nzt)[0]
    #elif mode == 'quad':
    #    sig = dblquad(lambda x, z: get_dsig(
    #        x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]
    return sig

#def get_sigP(xF, pT, rs, tar, had, mode='gauss', nx=100, nz=100):
#    xT = 2. * pT / rs
#    xF2 = xF * xF
#    xT2 = xT * xT
#
#    # Mandelstam variables at the hadron level
#    ss = rs * rs
#    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
#    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)
#
#    # Lower limits of the z and x integrations
#    zmin = np.sqrt(xF2 + xT2)
#
#    def xmin(z): return -uu / (z * ss + tt)
#
#    if mode == 'gauss':
#        dsigdzdx = np.vectorize(
#            lambda x, z: get_dsigP( x, z, xF, pT, rs, tar, had))
#        dsigdz = np.vectorize(lambda z: fixed_quad(
#            lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
#        sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
#    elif mode == 'quad':
#        sig = dblquad(lambda x, z: get_dsigP(
#            x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]
#    return sig

#def get_sigST(xF, pT, rs, tar, had, mode='gauss', nx=10, nz=10):
#    xT = 2. * pT / rs
#    xF2 = xF * xF
#    xT2 = xT * xT
#
#    # Mandelstam variables at the hadron level
#    ss = rs * rs
#    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
#    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)
#
#    # Lower limits of the z and x integrations
#    zmin = np.sqrt(xF2 + xT2)
#
#    def xmin(z): return -uu / (z * ss + tt)
#
#    if mode == 'gauss':
#      dsigdzdx = np.vectorize(
#        lambda x, z: get_dsigST(x, z, xF, pT, rs, tar, had))
#      dsigdz = np.vectorize(lambda z: fixed_quad(
#        lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
#      sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
#    elif mode == 'quad':
#      sig = dblquad(lambda x, z: get_dsigST(
#        x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]
#    return sig


if __name__ == '__main__':

    from qcdlib.lff0 import FF as FF0
    #from qcdlib.ff1 import FF as FF1
    from qcdlib.ldf0 import LDF as LDF0
    #from qcdlib.pdf1 import PDF as PDF1
    conf['aux']= AUX()
    conf['ldf']=LDF0()
    ##conf['collinspi']=FF1('pi')
    #conf['collinsk']=FF1('k')
    #conf['Htildepi']=FF1('pi')
    #conf['Htildek']=FF1('k')
    #conf['transversity']=PDF1()
    conf['ffmu']=LFF0('mu')
    conf['ffe']=LFF0('e')
    #conf['sivers']=PDF1()

    rs = 200.
    tar = 'p'
    had = 'pi-'
    pT = 1.10
    lpT = 1.10
    ppT = 1.10
    yL = .9
    yP = .9
    thLP = .785 #pi/4

    #xF = 0.2375
    #xT = 2. * pT / rs
    #xF2 = xF * xF
    #xT2 = xT * xT

    # Mandelstam variables at the hadron level
    ss = rs * rs
    #tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    #uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    # Lower limits of the z and x integrations
    #zmin = np.sqrt(xF2 + xT2)

    #def xmin(z): return -uu / (z * ss + tt)

    def test():
        den = get_sig(lpT,ppT,yL,yP,thLP, pT, rs, tar, had, mode='gauss', nx=100, nz=100)
        # num = get_sigST(xF, pT, rs, tar, had,
        #                      mode='gauss', nx=100, nz=100)
        # num2 = get_sigP(xF, pT, rs, tar, had,
        #                      mode='gauss', nx=100, nz=100)

        # AN = (num + num2) / den
        print den

#  from timeit import Timer
#  t = Timer("test()", "from __main__ import test")
#  print 't elapsed ',t.timeit(number=1)

#  def test2():
#    den = anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
#    num = anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
#
#    print den,num
#
#  from timeit import Timer
#  t = Timer("test2()", "from __main__ import test2")
#  print 't elapsed ',t.timeit(number=1)

#  start = time.time()
#  print anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
#  print anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
#  end = time.time()
#  print 'time=',(end-start)

    start = time.time()
    test()
    end = time.time()
    print 'time=', (end - start)

    # Integration of the numerator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)

    # Integration of the denominator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)
    #den = dblquad(lambda x,z: ANTHEORY().get_dsig(x,z,xF,pT,rs,tar,had),zmin,1.,xmin,lambda x: 1.)

    #AN = num[0]/den[0]
    # print(AN)
