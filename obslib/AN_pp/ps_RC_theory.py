#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import math
import time
from tools.tools import load
from qcdlib.aux import AUX
from qcdlib.alphaS import ALPHAS
from tools.config import conf, load_config
from scipy.integrate import quad, dblquad, fixed_quad

#ps_RC: calculate parton scattering unpolarized cross-section
class PSRC:
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


    def __init__(self,tag):

        load_config('data/input-%s.py'%tag)


        if 'basis' not in conf:
            conf['basis'] = 'default'

    def get_f(self,x, Q2): # Collinear unpolarized PDF
      return conf['pdf'].get_C(x, Q2)

    #def get_ft(self,x, Q2): # Collinear unpolarized PDF
    #  return conf['pdf'].get_C(x, Q2)

    def get_d(self,z, Q2, had): # Collinear unpolarized FF
      if 'pi' in had:
          return conf['ffpi'].get_C(z, Q2)
      elif 'k' in had:
          return conf['ffk'].get_C(z, Q2)

    def get_h(self,x, Q2): # Collinear transversity
      return conf['transversity'].get_C(x, Q2)

    def get_fl(self,xi, Q2): #LPDF
        #return 1.
        return conf['ldf'].get_C(xi, Q2)

    def get_dl(self,zt, Q2): #LFF
        #return 1.
        return conf['lff'].get_C(zt, Q2)

    def get_mandelstam(self,s, t, u, sp, tp, up):
    # Convenient combinations of the partonic Mandelstam variables
       m = {}
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

    def get_rapidityf(self,rs,ppT,lpT,yL,yP,thLP):
        rf = {}
        rf['ttp']=-rs*ppT*math.exp(yP)
        rf['uup']=-rs*lpT*math.exp(yL)
        rf['tt']=rs*lpT*math.exp(-yL)
        rf['uu']=-rs*ppT*math.exp(-yP)
        rf['ss']=rs*rs
        rf['ssp']=ppT*lpT*(math.exp(yL-yP)+math.exp(yP-yL)-2*math.cos(thLP))
        return rf

    def get_mparton(self,x,xi,z,zt,rs,ppT,lpT,yL,yP,thLP):
        mp={}

        oz = 1. / z
        ozt= 1. / zt

        rf=self.get_rapidityf(rs,ppT,lpT,yL,yP,thLP)

        mp['s'] = x * xi * rf['ss']
        mp['t'] = xi*ozt*rf['tt']
        mp['u'] = x * rf['uu'] * ozt
        mp['sp']  = oz*ozt*rf['ssp']
        mp['tp']  = x*oz*rf['ttp']
        mp['up']  = xi*oz*rf['uup']

        mp['Q2']=-rf['tt']

        return mp

    def get_Hupol(self,ppT,lpT,yL,yP,thLP,x, z, xi, zt, rs):
      # Hard parts for the unpolarized cross section

      #rf=get_rapidityf(rs,ppT,lpT,yL,yP,thLP)
      Hupol = {}

      mp=self.get_mparton(x,xi,z,zt,rs,ppT,lpT,yL,yP,thLP)

      m=self.get_mandelstam(mp['s'], mp['t'], mp['u'], mp['sp'], mp['tp'], mp['up'])


      Hupol[1] = 2*mp['Q2']*mp['Q2']*z*z
      Hupol[2] = 2*zt*mp['Q2']*mp['up']*z
      Hupol[3] = zt*zt*(m['s2']*x*x*z*z+2*m['up2'])
      Hupol[4] = 2*xi*mp['u']*x*z*(mp['Q2']*z-zt*mp['up'])
      Hupol[5] = m['u2']*x*x*z*z
      return Hupol


    # Calculation of the unpolarized cross section
    def get_dsig(self,ppT,lpT,yL,yP,thLP,x, z, xi, zt, rs, tar, had):

      # Mandelstam variables at the hadron level
      ss = rs * rs

      oz = 1. / z
      ozt= 1. / zt

      rf=self.get_rapidityf(rs,ppT,lpT,yL,yP,thLP)

      mp=self.get_mparton(x,xi,z,zt,rs,ppT,lpT,yL,yP,thLP)

      m=self.get_mandelstam(mp['s'], mp['t'], mp['u'], mp['sp'], mp['tp'], mp['up'])

      # Prefactor
      aS=ALPHAS()
      denfacint = 1. / (zt*zt*xi*z*z*x)
      pifac = (1. / math.pi**3.)
      denfac = 1. / (-xi*xi*rf['tt']*rf['ss']*rf['ss']*rf['uu']*x*x*z*z)
      numfac=-(1./137.)*as.get_alphaS(mp['Q2'])
      convert=(.389379)**3. #Gev^-6 to millibarn^3
      prefac=denfacint*pifac*denfac*numfac*convert

      Hupol=self.get_Hupol(ppT,lpT,yL,yP,thLP,x, z, xi, zt, rs)

      Hupol1 = Hupol[1]
      Hupol2 = Hupol[2]
      Hupol3 = Hupol[3]
      Hupol4 = Hupol[4]
      Hupol5 = Hupol[5]

      # Get arrays of the nonperturbative functions
      f = self.get_f(x, mp['Q2'])
      #ft = get_ft(xp, Q2)
      d = self.get_d(z, mp['Q2'], 'pi+')

      fl = self.get_fl(xi,mp['Q2'])

      dl = self.get_dl(zt,mp['Q2'])

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

      upol=fl*dl*Hupol*Hadprod

      return prefac * upol



    def get_sig(self,lpT,ppT,yL,yP,thLP, rs, tar, had,nz,nxi,nzt):


        rf=self.get_rapidityf(rs,ppT,lpT,yL,yP,thLP)

        # Lower limits of the  integrations
        ztmin = (-rf['tt']-rf['ssp']-rf['uup'])/(rf['ttp']+rf['ss']+rf['uu'])

        def ximin(zt): return (rf['ttp']+(1./zt)*rf['ssp']+(1./zt)*rf['uup'])/((1./zt)*rf['tt']-rf['ss']-rf['uu'])

        def zmin(xi,zt): return (rf['ttp']+(1/zt)*rf['ssp']+xi*rf['uu'])/((xi/zt)*rf['tt']-(1/zt)*rf['uup']-xi*rf['ss'])
    #x fixed from delta function??
        def x(xi,z,zt): return ((xi/zt)*rf['tt']-(1/(z*zt))*rf['ssp']+(xi/z)*rf['uu'])/((1/zt)*rf['uup']+xi*rf['ss']+(1/z)*rf['ttp'])


        dsigdztdxidz = np.vectorize(
            lambda z, xi, zt: self.get_dsig(ppT,lpT,yL,yP,thLP,x(xi,z,zt), z, xi,zt, rs, tar, had))
        dsigdztdxi = np.vectorize(lambda xi, zt: fixed_quad(
            lambda z: dsigdztdxidz(z, xi, zt), zmin(xi,zt), 1, n=nz)[0])
        dsigdzt = np.vectorize(lambda zt: fixed_quad(
            lambda xi: dsigdztdxi(xi, zt), ximin(zt), 1, n=nxi)[0])
        sig = fixed_quad(dsigdzt, ztmin, 1, n=nzt)[0]

        return sig


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
    #pT = 1.10
    lpT = 1.10
    ppT = 1.10
    yL = .9
    yP = .9
    thLP = .785 #pi/4


    ss = rs * rs


    def test():
        den = get_sig(lpT,ppT,yL,yP,thLP,rs, tar, had,nz=100,nxi=100, nzt=100)

        print den



    start = time.time()
    test()
    end = time.time()
    print 'time=', (end - start)
