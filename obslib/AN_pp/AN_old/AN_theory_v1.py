#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import math
import re
import time
from tools.tools import load_config
from external.CJLIB.CJ import CJ
from external.LSSLIB.LSS import LSS
from external.DSSLIB.DSS import DSS
from qcdlib.tmdlib import PDF, PPDF, FF
from qcdlib.tmdlib import TRANSVERSITY
from qcdlib.tmdlib import COLLINS
from qcdlib.aux import AUX
from tools.config import conf
from scipy.integrate import quad, dblquad, fixed_quad


class ANTHEORY:
    """
    AN_theory.py - program to calculate A_N in pp -> hX
    This only includes the fragmentation term (see 1701.09170)
    """
#    @profile

    def __init__(self):
        # self.aux=conf['aux']

        self.Mh = {}
        self.Mh['pi+'] = 0.13957018  # self.aux.Mpi
        self.Mh['pi-'] = 0.13957018  # self.aux.Mpi
        self.Mh['pi0'] = 0.1349766  # self.aux.Mpi
        self.Mh['k+'] = 0.13957018  # self.aux.Mk
        self.Mh['k-'] = 0.13957018  # self.aux.Mk

        self.flavor = ['u', 'd', 's', 'ub', 'db', 'sb', 'g']
        self.target = ['p']
        self.hadron = ['pi+', 'pi-', 'pi0']

        self.flavdict = {'g': 0, 'u': 1, 'd': 2,
                         's': 3, 'ub': 4, 'db': 5, 'sb': 6}

        # Common color factors and fractions
        self.c = {'r3': 1. / 3., 'r4': 0.25, 'r6': 1. / 6., 'r8': 0.125,
                  'r9': 1. / 9., 'r18': 1. / 18., 'r24': 1. / 24., 'r27': 1. / 27.}

        self.ch = {}
        self.chT = {}
        self.m = {}
        self.Hupol = {}
        self.f = {}
        self.ft = {}
        self.d = {}
        self.h = {}
        self.H1p = {}
        self.H = {}
        self.HTffa = np.zeros(12)
        self.HTffb = np.zeros(12)
        self.Hxxpz = np.zeros((12, 7))

        # Unpolarized cross section as strings to create matrices for array multiplication of pdf*pdf*ff*H
        #D=D1z, Ft=f1xp, F=f1x
        self.ch[0] = '(DdFdFtd+DsFsFtd+DdFdFts+DsFsFts+DdFdFtu+DsFsFtu+DuFuFtd+DuFuFts+DuFuFtu)*av1'
        self.ch[1] = '(DdFdFtd+DdFsFtd+DsFdFts+DsFsFts+DuFdFtu+DuFsFtu+DdFuFtd+DsFuFts+DuFuFtu)*av2'
        self.ch[2] = '(DdFdFtd+DsFsFts+DuFuFtu)*av3'
        self.ch[3] = '(DdbFdbFtdb+DsbFsbFtdb+DdbFdbFtsb+DsbFsbFtsb+DdbFdbFtub+DsbFsbFtub+DubFubFtdb+DubFubFtsb+DubFubFtub)*av1'
        self.ch[4] = '(DdbFdbFtdb+DdbFsbFtdb+DsbFdbFtsb+DsbFsbFtsb+DubFdbFtub+DubFsbFtub+DdbFubFtdb+DsbFubFtsb+DubFubFtub)*av2'
        self.ch[5] = '(DdbFdbFtdb+DsbFsbFtsb+DubFubFtub)*av3'
        self.ch[6] = '(DdFdFtdb+DsFsFtdb+DdFdFtsb+DsFsFtsb+DdFdFtub+DsFsFtub+DuFuFtdb+DuFuFtsb+DuFuFtub)*av5'
        self.ch[7] = '(DdFdFtdb+DsFdFtdb+DuFdFtdb+DdFsFtsb+DsFsFtsb+DuFsFtsb+DdFuFtub+DsFuFtub+DuFuFtub)*av4'
        self.ch[8] = '(DdFdFtdb+DsFsFtsb+DuFuFtub)*av6'
        self.ch[9] = '(DdbFdbFtd+DsbFsbFtd+DdbFdbFts+DsbFsbFts+DdbFdbFtu+DsbFsbFtu+DubFubFtd+DubFubFts+DubFubFtu)*av5'
        self.ch[10] = '(DdbFdbFtd+DsbFdbFtd+DubFdbFtd+DdbFsbFts+DsbFsbFts+DubFsbFts+DdbFubFtu+DsbFubFtu+DubFubFtu)*av4'
        self.ch[11] = '(DdbFdbFtd+DsbFsbFts+DubFubFtu)*av6'
        self.ch[12] = '(DdbFdFtdb+DdbFsFtdb+DsbFdFtsb+DsbFsFtsb+DubFdFtub+DubFsFtub+DdbFuFtdb+DsbFuFtsb+DubFuFtub)*av8'
        self.ch[13] = '(DdbFdFtdb+DsbFdFtdb+DubFdFtdb+DdbFsFtsb+DsbFsFtsb+DubFsFtsb+DdbFuFtub+DsbFuFtub+DubFuFtub)*av7'
        self.ch[14] = '(DdbFdFtdb+DsbFsFtsb+DubFuFtub)*av9'
        self.ch[15] = '(DdFdbFtd+DdFsbFtd+DsFdbFts+DsFsbFts+DuFdbFtu+DuFsbFtu+DdFubFtd+DsFubFts+DuFubFtu)*av8'
        self.ch[16] = '(DdFdbFtd+DsFdbFtd+DuFdbFtd+DdFsbFts+DsFsbFts+DuFsbFts+DdFubFtu+DsFubFtu+DuFubFtu)*av7'
        self.ch[17] = '(DdFdbFtd+DsFsbFts+DuFubFtu)*av9'
        self.ch[18] = '(DgFdFtdb+DgFsFtsb+DgFuFtub)*av10'
        self.ch[19] = '(DgFdbFtd+DgFsbFts+DgFubFtu)*av10'
        self.ch[20] = '(DdFdFtg+DsFsFtg+DuFuFtg)*av11'
        self.ch[21] = '(DdbFdbFtg+DsbFsbFtg+DubFubFtg)*av11'
        self.ch[22] = '(DgFdFtg+DgFsFtg+DgFuFtg)*av12'
        self.ch[23] = '(DgFdbFtg+DgFsbFtg+DgFubFtg)*av12'
        self.ch[24] = '(DdFgFtg+DsFgFtg+DuFgFtg)*av13'
        self.ch[25] = '(DdbFgFtg+DsbFgFtg+DubFgFtg)*av13'
        self.ch[26] = '(DgFgFtg)*av14'
        self.ch[27] = '(DdFgFtd+DsFgFts+DuFgFtu)*av12'
        self.ch[28] = '(DdbFgFtdb+DsbFgFtsb+DubFgFtub)*av12'
        self.ch[29] = '(DgFgFtd+DgFgFts+DgFgFtu)*av11'
        self.ch[30] = '(DgFgFtdb+DgFgFtsb+DgFgFtub)*av11'

        # Transversely polarized cross section as strings to create matrices for array multiplication of pdf*h1*Hxxpz
        #Ft=f1xp, HH=Hxxpz, h=h1x
        self.chT[0] = 'FtgHH0dhd+FtgHH0shs+FtgHH0uhu'
        self.chT[1] = 'FtdHH1dhd+FtsHH1shs+FtuHH1uhu'
        self.chT[2] = 'FtdbHH2dhd+FtsbHH2shs+FtubHH2uhu'
        self.chT[3] = 'FtdHH3dhdb+FtsHH3shsb+FtuHH3uhub'
        self.chT[4] = 'FtsHH4dhd+FtuHH4dhd+FtdHH4shs+FtuHH4shs+FtdHH4uhu+FtsHH4uhu'
        self.chT[5] = 'FtsbHH5dhd+FtubHH5dhd+FtdbHH5shs+FtubHH5shs+FtdbHH5uhu+FtsbHH5uhu'
        self.chT[6] = 'FtgHH6dbhdb+FtgHH6sbhsb+FtgHH6ubhub'
        self.chT[7] = 'FtdbHH7dbhdb+FtsbHH7sbhsb+FtubHH7ubhub'
        self.chT[8] = 'FtdHH8dbhdb+FtsHH8sbhsb+FtuHH8ubhub'
        self.chT[9] = 'FtdbHH9dbhd+FtsbHH9sbhs+FtubHH9ubhu'
        self.chT[10] = 'FtsbHH10dbhdb+FtubHH10dbhdb+FtdbHH10sbhsb+FtubHH10sbhsb+FtdbHH10ubhub+FtsbHH10ubhub'
        self.chT[11] = 'FtsHH11dbhdb+FtuHH11dbhdb+FtdHH11sbhsb+FtuHH11sbhsb+FtdHH11ubhub+FtsHH11ubhub'

        self.B = []
        for j in range(len(self.ch)):
            M = np.zeros((7, 7, 7))
            cs = self.ch[j].split('(')[1].split(')')[0].split('+')
            for i in range(len(cs)):
                _ = re.split('[DFFt]', cs[i])
                _ = filter(None, _)
                M[self.flavdict[_[0]], self.flavdict[_[1]], self.flavdict[_[2]]] = 1
            self.B.append(M)

        self.BT = []
        for j in range(len(self.chT)):
            MT = np.zeros((7, 12, 7, 7))
            csT = self.chT[j].split('+')
            for i in range(len(csT)):
                _ = re.split('[FtHHh]', csT[i])
                _ = filter(None, _)
                HHind = re.split('(\d+)', _[1])
                HHind = filter(None, HHind)
                MT[self.flavdict[_[0]], int(
                    HHind[0]), self.flavdict[HHind[1]], self.flavdict[_[2]]] = 1
            self.BT.append(MT)

    def get_f(self, x, Q2):  # Collinear unpolarized PDF
        if (x, Q2) not in self.f:
            self.f[(x, Q2)] = x * (1. - x) * np.ones(7)
        return self.f[(x, Q2)]

    def get_ft(self, x, Q2):  # Collinear unpolarized PDF
        if (x, Q2) not in self.ft:
            self.ft[(x, Q2)] = x * (1. - x) * np.ones(7)
        return self.ft[(x, Q2)]

    def get_d(self, z, Q2):  # Collinear unpolarized FF
        if (z, Q2) not in self.d:
            self.d[(z, Q2)] = z * (1. - z) * np.ones(7)
        return self.d[(z, Q2)]

    def get_h(self, x, Q2):  # Collinear transversity
        if (x, Q2) not in self.h:
            self.h[(x, Q2)] = x * (1. - x) * np.ones(7)
        return self.h[(x, Q2)]

    def get_H1p(self, z, Q2):  # (H_1^{\perp(1)}(z) - z*dH_1^{\perp(1)}(z)/dz)
        if (z, Q2) not in self.H1p:
            self.H1p[(z, Q2)] = z * (1. - z) * np.ones(7)
        return self.H1p[(z, Q2)]

    def get_H(self, z, Q2):  # -2*z*H_1^{\perp(1)}(z)+\tilde{H}(z)
        if (z, Q2) not in self.H:
            self.H[(z, Q2)] = z * (1. - z) * np.ones(7)
        return self.H[(z, Q2)]

    def get_mandelstam(self, s, t, u):
        # Convenient combinations of the partonic Mandelstam variables
        self.m['s2'] = s * s
        self.m['s3'] = s**3.
        self.m['t2'] = t * t
        self.m['t3'] = t**3.
        self.m['u2'] = u * u
        self.m['u3'] = u**3.
        self.m['ostu'] = 1. / (s * t * u)
        self.m['os'] = 1. / s
        self.m['ot'] = 1. / t
        self.m['ou'] = 1. / u
        self.m['st'] = s / t
        self.m['su'] = s / u
        self.m['ts'] = t / s
        self.m['tu'] = t / u
        self.m['us'] = u / s
        self.m['ut'] = u / t
        self.m['st2'] = s**2. / t**2.
        self.m['su2'] = s**2. / u**2.
        self.m['ts2'] = t**2. / s**2.
        self.m['tu2'] = t**2. / u**2.
        self.m['us2'] = u**2. / s**2.
        self.m['ut2'] = u**2. / t**2.
        self.m['os2'] = 1. / s**2.
        self.m['ot2'] = 1. / t**2.
        self.m['ou2'] = 1. / u**2.
        self.m['os3'] = 1. / s**3.
        self.m['ot3'] = 1. / t**3.
        self.m['ou3'] = 1. / u**3.

    def get_Hupol(self):
        # Hard parts for the unpolarized cross section
        self.Hupol['av1'] = 4. * (self.m['st2'] + self.m['ut2']) * self.c['r9']
        self.Hupol['av2'] = 4. * (self.m['su2'] + self.m['tu2']) * self.c['r9']
        self.Hupol['av3'] = -8. * self.m['st'] * self.m['su'] * self.c['r27']
        self.Hupol['av4'] = 4. * (self.m['ts2'] + self.m['us2']) * self.c['r9']
        self.Hupol['av5'] = 4. * (self.m['st2'] + self.m['ut2']) * self.c['r9']
        self.Hupol['av6'] = -8. * self.m['ut'] * self.m['us'] * self.c['r27']
        self.Hupol['av7'] = 4. * (self.m['ts2'] + self.m['us2']) * self.c['r9']
        self.Hupol['av8'] = 4. * (self.m['su2'] + self.m['tu2']) * self.c['r9']
        self.Hupol['av9'] = -8. * self.m['tu'] * self.m['ts'] * self.c['r27']
        self.Hupol['av10'] = 32. * self.c['r27'] * \
            (self.m['tu'] + self.m['ut']) * \
            (1. - 9. * self.m['ts'] * self.m['us'] * self.c['r4'])
        self.Hupol['av11'] = 4. * self.c['r9'] * (-self.m['su'] - self.m['us']) * (
            1. - 9. * self.m['st'] * self.m['ut'] * self.c['r4'])
        self.Hupol['av12'] = 4. * self.c['r9'] * \
            (-self.m['st'] - self.m['ts']) * \
            (1. - 9. * self.m['su'] * self.m['tu'] * self.c['r4'])
        self.Hupol['av13'] = (self.m['tu'] + self.m['ut']) * self.c['r6'] - \
            3. * (self.m['ts2'] + self.m['us2']) * self.c['r8']
        self.Hupol['av14'] = 4.5 * (3. - self.m['ts'] * self.m['us'] -
                                  self.m['st'] * self.m['ut'] - self.m['su'] * self.m['tu'])

    def get_HTffa(self, s, t, u):
        # Hard parts for the transversely polarized fragmentation term
        self.HTffa[0] = -self.c['r9'] * self.m['ot'] + self.c['r8'] * \
            s * (u - s) * self.m['ot3'] - self.m['st2'] * self.m['ou']

        self.HTffa[1] = self.c['r27'] * s * (t - u) * self.m['ot2'] * self.m['ou'] + \
            self.c['r9'] * s * (u - 2. * t) * self.m['ot3'] + s * self.m['ot2']

        self.HTffa[2] = self.c['r27'] * s * self.m['ot2'] + self.c['r9'] * \
            s * (t - s) * self.m['ot3'] - self.c['r3'] * self.m['ot']

        self.HTffa[3] = self.c['r27'] * s * self.m['ot'] * \
            self.m['ou'] - self.c['r3'] * self.m['ot']

        self.HTffa[4] = self.c['r9'] * s * (u - 2. * t) * self.m['ot3'] + s * self.m['ot2']

        self.HTffa[5] = self.c['r9'] * s * (t - s) * self.m['ot3']

        self.HTffa[6] = self.HTffa[0]

        self.HTffa[7] = self.HTffa[1]

        self.HTffa[8] = self.HTffa[2]

        self.HTffa[9] = self.HTffa[3]

        self.HTffa[10] = self.HTffa[4]

        self.HTffa[11] = self.HTffa[5]

        return self.HTffa

    def get_HTffb(self, s, t, u):
        # Hard parts for the transversely polarized fragmentation term
        self.HTffb[0] = self.c['r8'] * s * (u - s) * self.m['ot3'] + 0.5 * self.c['r9'] * (s - u) * self.m['ot'] * \
            self.m['ou'] + 0.5 * (s - u) * (self.m['t2'] - 2. * t *
                                      u - 2. * self.m['u2']) * self.m['ot3'] * self.m['ou']

        self.HTffb[1] = self.c['r27'] * 0.5 * s * (t - 3. * u) * self.m['ot2'] * self.m['ou'] - s * u * self.m['ot3'] + self.c['r9'] * s * (
            2. * u - t) * self.m['ot3'] - self.c['r3'] * 0.5 * self.m['s2'] * self.m['ot2'] * self.m['ou']

        self.HTffb[2] = self.c['r27'] * 0.5 * (3. * s - t) * self.m['ot2'] + self.m['s2'] * self.m['ot3'] + \
            self.c['r9'] * s * (t - 2. * s) * self.m['ot3'] + \
            self.c['r3'] * 0.5 * u * self.m['ot2']

        self.HTffb[3] = 10. * self.c['r27'] * 0.5 * (s - u) * self.m['ot'] * self.m['ou']

        self.HTffb[4] = self.c['r9'] * s * \
            (2. * u - t) * self.m['ot3'] - s * u * self.m['ot3']

        self.HTffb[5] = self.c['r9'] * s * \
            (t - 2. * s) * self.m['ot3'] + self.m['s2'] * self.m['ot3']

        self.HTffb[6] = self.HTffb[0]

        self.HTffb[7] = self.HTffb[1]

        self.HTffb[8] = self.HTffb[2]

        self.HTffb[9] = self.HTffb[3]

        self.HTffb[10] = self.HTffb[4]

        self.HTffb[11] = self.HTffb[5]

        return self.HTffb

    def get_Hxxpz(self, z, Q2, s, t, u):

        HTffa = self.get_HTffa(s, t, u)
        HTffb = self.get_HTffb(s, t, u)
        H1p = self.get_H1p(z, Q2)
        H = self.get_H(z, Q2)

        self.Hxxpz = np.einsum('i,j->ij', HTffa, H1p) + \
            np.einsum('i,j->ij', HTffb, H) / z

        return self.Hxxpz

#    @profile
    # Calculation of the unpolarized cross section
    def get_dsig(self, x, z, xF, pT, rs, tar, had):

        if pT < 1.:
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

        xp = -x * tt / (z * x * ss + uu)

        # Mandelstam variables at the parton level
        s = x * xp * ss
        t = x * tt * oz
        u = xp * uu * oz

        # Prefactor
        denfac = 1. / ((z * z * x * ss + uu * z) * x * xp)

        self.get_mandelstam(s, t, u)
        self.get_Hupol()
#        self.get_matrices()

        # Get arrays of the nonperturbative functions
        f = self.get_f(x, Q2)
        ft = self.get_ft(xp, Q2)
        d = self.get_d(x, Q2)

        upol = 0
        for i in range(len(self.ch)):
            upol += self.Hupol[self.ch[i].split('*')[1]] * \
                np.einsum('i,j,k,ijk', d, f, ft, self.B[i])

        return denfac * upol

#    @profile
    # Calculation of the fragmentation term in the transversely polarized cross section
    def get_dsigST(self, x, z, xF, pT, rs, tar, had):

        Mh = self.Mh[had]

        if pT < 1.:
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

        xp = -x * tt / (z * x * ss + uu)

        # Mandelstam variables at the parton level
        s = x * xp * ss
        t = x * tt * oz
        u = xp * uu * oz

        # Prefactor
        numfac = oz * 1. / ((z * z * x * ss + uu * z) * x * xp)

        self.get_mandelstam(s, t, u)

        # Get arrays of the nonperturbative functions
        ft = self.get_ft(xp, Q2)
        h = self.get_h(x, Q2)
        Hxxpz = self.get_Hxxpz(z, Q2, s, t, u)

        ffcs = 0
        for i in range(len(self.chT)):
            ffcs += np.einsum('i,jk,l,ijkl', ft, Hxxpz, h, self.BT[i])

        ffcs = 2. * Mh * pT * numfac * ffcs

        return ffcs

    def get_sig(self, xF, pT, rs, tar, had, mode='gauss', nx=100, nz=100):
        if mode == 'gauss':
            dsigdzdx = np.vectorize(
                lambda x, z: self.get_dsig(x, z, xF, pT, rs, tar, had))
            dsigdz = np.vectorize(lambda z: fixed_quad(
                lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
            sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
        elif mode == 'quad':
            sig = dblquad(lambda x, z: self.get_dsig(
                x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]
        return sig

    def get_sigST(self, xF, pT, rs, tar, had, mode='gauss', nx=100, nz=100):
        if mode == 'gauss':
            dsigdzdx = np.vectorize(
                lambda x, z: self.get_dsigST(x, z, xF, pT, rs, tar, had))
            dsigdz = np.vectorize(lambda z: fixed_quad(
                lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
            sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
        elif mode == 'quad':
            sig = dblquad(lambda x, z: self.get_dsigST(
                x, z, xF, pT, rs, tar, had), zmin, 1., xmin, lambda x: 1.)[0]
        return sig


if __name__ == '__main__':

    rs = 200.
    tar = 'p'
    had = 'pi-'
    pT = 1.10
    xF = 0.2375
    xT = 2. * pT / rs
    xF2 = xF * xF
    xT2 = xT * xT

    # Mandelstam variables at the hadron level
    ss = rs * rs
    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    # Lower limits of the z and x integrations
    zmin = np.sqrt(xF2 + xT2)

    def xmin(z): return -uu / (z * ss + tt)

    anthy = ANTHEORY()

    def test():
        den = anthy.get_sig(xF, pT, rs, tar, had, mode='gauss', nx=100, nz=100)
        num = anthy.get_sigST(xF, pT, rs, tar, had,
                              mode='gauss', nx=100, nz=100)

        AN = num / den
        print AN
#
#  test()
#
#  from timeit import Timer
#  t = Timer("test()", "from __main__ import test")
#  print 't elapsed ',t.timeit(number=1)

#  start = time.time()
#  print anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
#  print anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
#  end = time.time()
#  print(end-start)

    start = time.time()
    test()
    end = time.time()
    print 'time=', (end - start)

    # Integration of the numerator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)

    # Integration of the denominator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)
    #den = dblquad(lambda x,z: ANTHEORY().get_dsig(x,z,xF,pT,rs,tar,had),zmin,1.,xmin,lambda x: 1.)

    #AN = num[0]/den[0]
    # print(AN)
