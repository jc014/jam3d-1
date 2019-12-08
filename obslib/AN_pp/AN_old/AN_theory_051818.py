#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import math
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


class ANTHEORY():
    """
    AN_theory.py - program to calculate A_N in pp -> hX
    This only includes the fragmentation term (see 1701.09170)
    """

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

        self.m = {}
        self.Hupol = {}
        self.f = {}
        self.ft = {}
        self.d = {}
        self.h = {}
        self.H1p = {}
        self.H = {}
        self.HTffa = np.zeros(13)
        self.HTffb = np.zeros(13)
        self.Hxxpz = np.zeros((13, 7))

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

        return self.m

    def get_Hupol(self):
        # Hard parts for the unpolarized cross section
        self.Hupol[1] = 4. * (self.m['st2'] + self.m['ut2']) * self.c['r9']
        self.Hupol[2] = 4. * (self.m['su2'] + self.m['tu2']) * self.c['r9']
        self.Hupol[3] = -8. * self.m['st'] * self.m['su'] * self.c['r27']
        self.Hupol[4] = 4. * (self.m['ts2'] + self.m['us2']) * self.c['r9']
        self.Hupol[5] = 4. * (self.m['st2'] + self.m['ut2']) * self.c['r9']
        self.Hupol[6] = -8. * self.m['ut'] * self.m['us'] * self.c['r27']
        self.Hupol[7] = 4. * (self.m['ts2'] + self.m['us2']) * self.c['r9']
        self.Hupol[8] = 4. * (self.m['su2'] + self.m['tu2']) * self.c['r9']
        self.Hupol[9] = -8. * self.m['tu'] * self.m['ts'] * self.c['r27']
        self.Hupol[10] = 32. * self.c['r27'] * \
            (self.m['tu'] + self.m['ut']) * \
            (1. - 9. * self.m['ts'] * self.m['us'] * self.c['r4'])
        self.Hupol[11] = 4. * self.c['r9'] * (-self.m['su'] - self.m['us']) * (
            1. - 9. * self.m['st'] * self.m['ut'] * self.c['r4'])
        self.Hupol[12] = 4. * self.c['r9'] * (-self.m['st'] - self.m['ts']) * (
            1. - 9. * self.m['su'] * self.m['tu'] * self.c['r4'])
        self.Hupol[13] = (self.m['tu'] + self.m['ut']) * self.c['r6'] - \
            3. * (self.m['ts2'] + self.m['us2']) * self.c['r8']
        self.Hupol[14] = 4.5 * (3. - self.m['ts'] * self.m['us'] -
                              self.m['st'] * self.m['ut'] - self.m['su'] * self.m['tu'])

        return self.Hupol

    def get_HTffa(self, s, t, u):
        # Hard parts for the transversely polarized fragmentation term
        self.HTffa[0] = 0

        self.HTffa[1] = -self.c['r9'] * self.m['ot'] + self.c['r8'] * \
            s * (u - s) * self.m['ot3'] - self.m['st2'] * self.m['ou']

        self.HTffa[2] = self.c['r27'] * s * (t - u) * self.m['ot2'] * self.m['ou'] + \
            self.c['r9'] * s * (u - 2. * t) * self.m['ot3'] + s * self.m['ot2']

        self.HTffa[3] = self.c['r27'] * s * self.m['ot2'] + self.c['r9'] * \
            s * (t - s) * self.m['ot3'] - self.c['r3'] * self.m['ot']

        self.HTffa[4] = self.c['r27'] * s * self.m['ot'] * \
            self.m['ou'] - self.c['r3'] * self.m['ot']

        self.HTffa[5] = self.c['r9'] * s * (u - 2. * t) * self.m['ot3'] + s * self.m['ot2']

        self.HTffa[6] = self.c['r9'] * s * (t - s) * self.m['ot3']

        self.HTffa[7] = self.HTffa[1]

        self.HTffa[8] = self.HTffa[2]

        self.HTffa[9] = self.HTffa[3]

        self.HTffa[10] = self.HTffa[4]

        self.HTffa[11] = self.HTffa[5]

        self.HTffa[12] = self.HTffa[6]

        return self.HTffa

    def get_HTffb(self, s, t, u):
        # Hard parts for the transversely polarized fragmentation term
        self.HTffb[0] = 0

        self.HTffb[1] = self.c['r8'] * s * (u - s) * self.m['ot3'] + 0.5 * self.c['r9'] * (s - u) * self.m['ot'] * \
            self.m['ou'] + 0.5 * (s - u) * (self.m['t2'] - 2. * t *
                                      u - 2. * self.m['u2']) * self.m['ot3'] * self.m['ou']

        self.HTffb[2] = self.c['r27'] * 0.5 * s * (t - 3. * u) * self.m['ot2'] * self.m['ou'] - s * u * self.m['ot3'] + self.c['r9'] * s * (
            2. * u - t) * self.m['ot3'] - self.c['r3'] * 0.5 * self.m['s2'] * self.m['ot2'] * self.m['ou']

        self.HTffb[3] = self.c['r27'] * 0.5 * (3. * s - t) * self.m['ot2'] + self.m['s2'] * self.m['ot3'] + \
            self.c['r9'] * s * (t - 2. * s) * self.m['ot3'] + \
            self.c['r3'] * 0.5 * u * self.m['ot2']

        self.HTffb[4] = 10. * self.c['r27'] * 0.5 * (s - u) * self.m['ot'] * self.m['ou']

        self.HTffb[5] = self.c['r9'] * s * \
            (2. * u - t) * self.m['ot3'] - s * u * self.m['ot3']

        self.HTffb[6] = self.c['r9'] * s * \
            (t - 2. * s) * self.m['ot3'] + self.m['s2'] * self.m['ot3']

        self.HTffb[7] = self.HTffb[1]

        self.HTffb[8] = self.HTffb[2]

        self.HTffb[9] = self.HTffb[3]

        self.HTffb[10] = self.HTffb[4]

        self.HTffb[11] = self.HTffb[5]

        self.HTffb[12] = self.HTffb[6]

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

        Hupol1 = self.Hupol[1]
        Hupol2 = self.Hupol[2]
        Hupol3 = self.Hupol[3]
        Hupol4 = self.Hupol[4]
        Hupol5 = self.Hupol[5]
        Hupol6 = self.Hupol[6]
        Hupol7 = self.Hupol[7]
        Hupol8 = self.Hupol[8]
        Hupol9 = self.Hupol[9]
        Hupol10 = self.Hupol[10]
        Hupol11 = self.Hupol[11]
        Hupol12 = self.Hupol[12]
        Hupol13 = self.Hupol[13]
        Hupol14 = self.Hupol[14]

        # Get arrays of the nonperturbative functions
        f = self.get_f(x, Q2)
        ft = self.get_ft(xp, Q2)
        d = self.get_d(x, Q2)

        fg = f[0]
        fu = f[1]
        fd = f[2]
        fs = f[3]
        fub = f[4]
        fdb = f[5]
        fsb = f[6]

        ftg = ft[0]
        ftu = ft[1]
        ftd = ft[2]
        fts = ft[3]
        ftub = ft[4]
        ftdb = ft[5]
        ftsb = ft[6]

        dg = d[0]
        du = d[1]
        dd = d[2]
        ds = d[3]
        dub = d[4]
        ddb = d[5]
        dsb = d[6]

        upol = 0

        upol += (fu * du + fd * dd + fs * ds) * (ftu + ftd + fts) * Hupol1
        + (fu + fd + fs) * (ftu * du + ftd * dd + fts * ds) * Hupol2
        + (fu * ftu * du + fd * ftd * dd + fs * fts * ds) * Hupol3

        upol += (fub * dub + fdb * ddb + fsb * dsb) * (ftub + ftdb + ftsb) * Hupol1
        + (fub + fdb + fsb) * (ftub * dub + ftdb * ddb + ftsb * dsb) * Hupol2
        + (fub * ftub * dub + fdb * ftdb * ddb + fsb * ftsb * dsb) * Hupol3

        upol += (fu * du + fd * dd + fs * ds) * (ftub + ftdb + ftsb) * Hupol5
        + (fu * ftub + fd * ftdb + fs * ftsb) * (du + dd + ds) * Hupol4
        + (fu * ftub * du + fd * ftdb * dd + fs * ftsb * ds) * Hupol6

        upol += (fub * dub + fdb * ddb + fsb * dsb) * (ftu + ftd + fts) * Hupol5
        + (fub * ftu + fdb * ftd + fsb * fts) * (dub + ddb + dsb) * Hupol4
        + (fub * ftu * dub + fdb * ftd * ddb + fsb * fts * dsb) * Hupol6

        upol += (fu + fd + fs) * (ftub * dub + ftdb * ddb + ftsb * dsb) * Hupol8
        + (fu * ftub + fd * ftdb + fs * ftsb) * (dub + ddb + dsb) * Hupol7
        + (fu * ftub * dub + fd * ftdb * ddb + fs * ftsb * dsb) * Hupol9

        upol += (fub + fdb + fsb) * (ftu * du + ftd * dd + fts * ds) * Hupol8
        + (fub * ftu + fdb * ftd + fsb * fts) * (du + dd + ds) * Hupol7
        + (fub * ftu * du + fdb * ftd * dd + fsb * fts * ds) * Hupol9

        upol += (fu * ftub + fd * ftdb + fs * ftsb) * dg * Hupol10

        upol += (fub * ftu + fdb * ftd + fsb * fts) * dg * Hupol10

        upol += (fu * du + fd * dd + fs * ds) * ftg * Hupol11

        upol += (fub * dub + fdb * ddb + fsb * dsb) * ftg * Hupol11

        upol += (fu + fd + fs) * ftg * dg * Hupol12

        upol += (fub + fdb + fsb) * ftg * dg * Hupol12

        upol += fg * ftg * (du + dd + ds) * Hupol13

        upol += fg * ftg * (dub + ddb + dsb) * Hupol13

        upol += fg * ftg * dg * Hupol14

        upol += fg * (ftu * du + ftd * dd + fts * ds) * Hupol12

        upol += fg * (ftub * dub + ftdb * ddb + ftsb * dsb) * Hupol12

        upol += fg * (ftu + ftd + fts) * dg * Hupol11

        upol += fg * (ftub + ftdb + ftsb) * dg * Hupol11

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

        hg = h[0]
        hu = h[1]
        hd = h[2]
        hs = h[3]
        hub = h[4]
        hdb = h[5]
        hsb = h[6]

        ftg = ft[0]
        ftu = ft[1]
        ftd = ft[2]
        fts = ft[3]
        ftub = ft[4]
        ftdb = ft[5]
        ftsb = ft[6]

        Hxxpz1u = Hxxpz[1][1]
        Hxxpz1d = Hxxpz[1][2]
        Hxxpz1s = Hxxpz[1][3]
        Hxxpz2u = Hxxpz[2][1]
        Hxxpz2d = Hxxpz[2][2]
        Hxxpz2s = Hxxpz[2][3]
        Hxxpz3u = Hxxpz[3][1]
        Hxxpz3d = Hxxpz[3][2]
        Hxxpz3s = Hxxpz[3][3]
        Hxxpz4u = Hxxpz[4][1]
        Hxxpz4d = Hxxpz[4][2]
        Hxxpz4s = Hxxpz[4][3]
        Hxxpz5u = Hxxpz[5][1]
        Hxxpz5d = Hxxpz[5][2]
        Hxxpz5s = Hxxpz[5][3]
        Hxxpz6u = Hxxpz[6][1]
        Hxxpz6d = Hxxpz[6][2]
        Hxxpz6s = Hxxpz[6][3]

        Hxxpz7ub = Hxxpz[7][1]
        Hxxpz7db = Hxxpz[7][2]
        Hxxpz7sb = Hxxpz[7][3]
        Hxxpz8ub = Hxxpz[8][1]
        Hxxpz8db = Hxxpz[8][2]
        Hxxpz8sb = Hxxpz[8][3]
        Hxxpz9ub = Hxxpz[9][1]
        Hxxpz9db = Hxxpz[9][2]
        Hxxpz9sb = Hxxpz[9][3]
        Hxxpz10ub = Hxxpz[10][1]
        Hxxpz10db = Hxxpz[10][2]
        Hxxpz10sb = Hxxpz[10][3]
        Hxxpz11ub = Hxxpz[11][1]
        Hxxpz11db = Hxxpz[11][2]
        Hxxpz11sb = Hxxpz[11][3]
        Hxxpz12ub = Hxxpz[12][1]
        Hxxpz12db = Hxxpz[12][2]
        Hxxpz12sb = Hxxpz[12][3]

        ffcs = 0

        ffcs += ftg * (hu * Hxxpz1u + hd * Hxxpz1d + hs * Hxxpz1s)

        ffcs += (hu * ftu * Hxxpz2u) + (hd * ftd * Hxxpz2d) + (hs * fts * Hxxpz2s)

        ffcs += (hu * ftub * Hxxpz3u) + (hd * ftdb * Hxxpz3d) + (hs * ftsb * Hxxpz3s)

        ffcs += (hub * ftu * Hxxpz4u) + (hdb * ftd * Hxxpz4d) + (hsb * fts * Hxxpz4s)

        ffcs += ftu * (hd * Hxxpz5d + hs * Hxxpz5s) + ftd * (hu *
                                                         Hxxpz5u + hs * Hxxpz5s) + fts * (hu * Hxxpz5u + hd * Hxxpz5d)

        ffcs += ftub * (hd * Hxxpz6d + hs * Hxxpz6s) + ftdb * (hu *
                                                           Hxxpz6u + hs * Hxxpz6s) + ftsb * (hu * Hxxpz6u + hd * Hxxpz6d)

        ffcs += ftg * (hub * Hxxpz7ub + hdb * Hxxpz7db + hsb * Hxxpz7sb)

        ffcs += (hub * ftub * Hxxpz8ub) + (hdb * ftdb * Hxxpz8db) + (hsb * ftsb * Hxxpz8sb)

        ffcs += (hub * ftu * Hxxpz9ub) + (hdb * ftd * Hxxpz9db) + (hsb * fts * Hxxpz9sb)

        ffcs += (hu * ftub * Hxxpz10ub) + (hd * ftdb * Hxxpz10db) + (hs * ftsb * Hxxpz10sb)

        ffcs += ftub * (hdb * Hxxpz11db + hsb * Hxxpz11sb) + ftdb * (hub *
                                                                 Hxxpz11ub + hsb * Hxxpz11sb) + fts * (hub * Hxxpz11ub + hdb * Hxxpz11db)

        ffcs += ftu * (hdb * Hxxpz12db + hsb * Hxxpz12sb) + ftd * (hub *
                                                               Hxxpz12ub + hsb * Hxxpz12sb) + fts * (hub * Hxxpz12ub + hdb * Hxxpz12db)

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
