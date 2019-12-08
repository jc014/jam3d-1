#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import math
import time
from numba import jit, float32
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


"""
AN_theory.py - program to calculate A_N in pp -> hX
This only includes the fragmentation term (see 1701.09170)
"""

flavor = ['u', 'd', 's', 'ub', 'db', 'sb', 'g']
target = ['p']
hadron = ['pi+', 'pi-', 'pi0']

flavdict = {'g': 0, 'u': 1, 'd': 2, 's': 3, 'ub': 4, 'db': 5, 'sb': 6}


@jit(float32(float32, float32))
def get_f(x, Q2):  # Collinear unpolarized PDF
    f = {}
    f[(x, Q2)] = x * (1. - x) * np.ones(7)
    return f[(x, Q2)]


@jit(float32(float32, float32))
def get_ft(x, Q2):  # Collinear unpolarized PDF
    ft = {}
    ft[(x, Q2)] = x * (1. - x) * np.ones(7)
    return ft[(x, Q2)]


@jit(float32(float32, float32))
def get_d(z, Q2):  # Collinear unpolarized FF
    d = {}
    d[(z, Q2)] = z * (1. - z) * np.ones(7)
    return d[(z, Q2)]


@jit(float32(float32, float32))
def get_h(x, Q2):  # Collinear transversity
    h = {}
    h[(x, Q2)] = x * (1. - x) * np.ones(7)
    return h[(x, Q2)]


@jit(float32(float32, float32))
def get_H1p(z, Q2):  # (H_1^{\perp(1)}(z) - z*dH_1^{\perp(1)}(z)/dz)
    H1p = {}
    H1p[(z, Q2)] = z * (1. - z) * np.ones(7)
    return H1p[(z, Q2)]


@jit(float32(float32, float32))
def get_H(z, Q2):  # -2*z*H_1^{\perp(1)}(z)+\tilde{H}(z)
    H = {}
    H[(z, Q2)] = z * (1. - z) * np.ones(7)
    return H[(z, Q2)]


@jit(float32[:](float32, float32, float32))
def get_Hupol(s, t, u):

    Hupol = {}

    # Common color factors and fractions
    r3 = 1. / 3.
    r4 = 0.25
    r6 = 1. / 6.
    r8 = 0.125
    r9 = 1. / 9.
    r18 = 1. / 18.
    r24 = 1. / 24.
    r27 = 1. / 27.

    # Convenient combinations of the partonic Mandelstam variables
    s2 = s * s
    s3 = s**3.
    t2 = t * t
    t3 = t**3.
    u2 = u * u
    u3 = u**3.
    ostu = 1. / (s * t * u)
    os = 1. / s
    ot = 1. / t
    ou = 1. / u
    st = s / t
    su = s / u
    ts = t / s
    tu = t / u
    us = u / s
    ut = u / t
    st2 = s**2. / t**2.
    su2 = s**2. / u**2.
    ts2 = t**2. / s**2.
    tu2 = t**2. / u**2.
    us2 = u**2. / s**2.
    ut2 = u**2. / t**2.
    os2 = 1. / s**2.
    ot2 = 1. / t**2.
    ou2 = 1. / u**2.
    os3 = 1. / s**3.
    ot3 = 1. / t**3.
    ou3 = 1. / u**3.

    # Hard parts for the unpolarized cross section
    Hupol[1] = 4. * (st2 + ut2) * r9
    Hupol[2] = 4. * (su2 + tu2) * r9
    Hupol[3] = -8. * st * su * r27
    Hupol[4] = 4. * (ts2 + us2) * r9
    Hupol[5] = 4. * (st2 + ut2) * r9
    Hupol[6] = -8. * ut * us * r27
    Hupol[7] = 4. * (ts2 + us2) * r9
    Hupol[8] = 4. * (su2 + tu2) * r9
    Hupol[9] = -8. * tu * ts * r27
    Hupol[10] = 32. * r27 * (tu + ut) * (1. - 9. * ts * us * r4)
    Hupol[11] = 4. * r9 * (-su - us) * (1. - 9. * st * ut * r4)
    Hupol[12] = 4. * r9 * (-st - ts) * (1. - 9. * su * tu * r4)
    Hupol[13] = (tu + ut) * r6 - 3. * (ts2 + us2) * r8
    Hupol[14] = 4.5 * (3. - ts * us - st * ut - su * tu)

    return Hupol


@jit(float32[:](float32, float32, float32))
def get_HTffa(s, t, u):

    HTffa = np.zeros(13)

    # Common color factors and fractions
    r3 = 1. / 3.
    r4 = 0.25
    r6 = 1. / 6.
    r8 = 0.125
    r9 = 1. / 9.
    r18 = 1. / 18.
    r24 = 1. / 24.
    r27 = 1. / 27.

    # Convenient combinations of the partonic Mandelstam variables
    s2 = s * s
    s3 = s**3.
    t2 = t * t
    t3 = t**3.
    u2 = u * u
    u3 = u**3.
    ostu = 1. / (s * t * u)
    os = 1. / s
    ot = 1. / t
    ou = 1. / u
    st = s / t
    su = s / u
    ts = t / s
    tu = t / u
    us = u / s
    ut = u / t
    st2 = s**2. / t**2.
    su2 = s**2. / u**2.
    ts2 = t**2. / s**2.
    tu2 = t**2. / u**2.
    us2 = u**2. / s**2.
    ut2 = u**2. / t**2.
    os2 = 1. / s**2.
    ot2 = 1. / t**2.
    ou2 = 1. / u**2.
    os3 = 1. / s**3.
    ot3 = 1. / t**3.
    ou3 = 1. / u**3.

    # Hard parts for the transversely polarized fragmentation term
    HTffa[0] = 0

    HTffa[1] = -r9 * ot + r8 * s * (u - s) * ot3 - st2 * ou

    HTffa[2] = r27 * s * (t - u) * ot2 * ou + r9 * s * (u - 2. * t) * ot3 + s * ot2

    HTffa[3] = r27 * s * ot2 + r9 * s * (t - s) * ot3 - r3 * ot

    HTffa[4] = r27 * s * ot * ou - r3 * ot

    HTffa[5] = r9 * s * (u - 2. * t) * ot3 + s * ot2

    HTffa[6] = r9 * s * (t - s) * ot3

    HTffa[7] = HTffa[1]

    HTffa[8] = HTffa[2]

    HTffa[9] = HTffa[3]

    HTffa[10] = HTffa[4]

    HTffa[11] = HTffa[5]

    HTffa[12] = HTffa[6]

    return HTffa


@jit(float32[:](float32, float32, float32))
def get_HTffb(s, t, u):

    HTffb = np.zeros(13)

    # Common color factors and fractions
    r3 = 1. / 3.
    r4 = 0.25
    r6 = 1. / 6.
    r8 = 0.125
    r9 = 1. / 9.
    r18 = 1. / 18.
    r24 = 1. / 24.
    r27 = 1. / 27.

    # Convenient combinations of the partonic Mandelstam variables
    s2 = s * s
    s3 = s**3.
    t2 = t * t
    t3 = t**3.
    u2 = u * u
    u3 = u**3.
    ostu = 1. / (s * t * u)
    os = 1. / s
    ot = 1. / t
    ou = 1. / u
    st = s / t
    su = s / u
    ts = t / s
    tu = t / u
    us = u / s
    ut = u / t
    st2 = s**2. / t**2.
    su2 = s**2. / u**2.
    ts2 = t**2. / s**2.
    tu2 = t**2. / u**2.
    us2 = u**2. / s**2.
    ut2 = u**2. / t**2.
    os2 = 1. / s**2.
    ot2 = 1. / t**2.
    ou2 = 1. / u**2.
    os3 = 1. / s**3.
    ot3 = 1. / t**3.
    ou3 = 1. / u**3.

    # Hard parts for the transversely polarized fragmentation term
    HTffb[0] = 0

    HTffb[1] = r8 * s * (u - s) * ot3 + 0.5 * r9 * (s - u) * ot * ou + \
        0.5 * (s - u) * (t2 - 2. * t * u - 2. * u2) * ot3 * ou

    HTffb[2] = r27 * 0.5 * s * (t - 3. * u) * ot2 * ou - s * u * ot3 + r9 * \
        s * (2. * u - t) * ot3 - r3 * 0.5 * s2 * ot2 * ou

    HTffb[3] = r27 * 0.5 * (3. * s - t) * ot2 + s2 * ot3 + r9 * s * (t - 2. * s) * ot3 + r3 * 0.5 * u * ot2

    HTffb[4] = 10. * r27 * 0.5 * (s - u) * ot * ou

    HTffb[5] = r9 * s * (2. * u - t) * ot3 - s * u * ot3

    HTffb[6] = r9 * s * (t - 2. * s) * ot3 + s2 * ot3

    HTffb[7] = HTffb[1]

    HTffb[8] = HTffb[2]

    HTffb[9] = HTffb[3]

    HTffb[10] = HTffb[4]

    HTffb[11] = HTffb[5]

    HTffb[12] = HTffb[6]

    return HTffb


def get_Hxxpz(z, Q2, s, t, u):

    Hxxpz = np.zeros((13, 7))

    HTffa = get_HTffa(s, t, u)
    HTffb = get_HTffb(s, t, u)
    H1p = get_H1p(z, Q2)
    H = get_H(z, Q2)

    Hxxpz = np.einsum('i,j->ij', HTffa, H1p) + np.einsum('i,j->ij', HTffb, H) / z

    return Hxxpz

#   @profile


@jit(float32(float32, float32, float32, float32, float32))
def get_dsig(x, z, xF, pT, rs):  # Calculation of the unpolarized cross section

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

    Hupol = get_Hupol(s, t, u)

    Hupol1 = Hupol[1]
    Hupol2 = Hupol[2]
    Hupol3 = Hupol[3]
    Hupol4 = Hupol[4]
    Hupol5 = Hupol[5]
    Hupol6 = Hupol[6]
    Hupol7 = Hupol[7]
    Hupol8 = Hupol[8]
    Hupol9 = Hupol[9]
    Hupol10 = Hupol[10]
    Hupol11 = Hupol[11]
    Hupol12 = Hupol[12]
    Hupol13 = Hupol[13]
    Hupol14 = Hupol[14]

    # Get arrays of the nonperturbative functions
    f = get_f(x, Q2)
    ft = get_ft(xp, Q2)
    d = get_d(x, Q2)

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

#   @profile


@jit(float32(float32, float32, float32, float32, float32))
# Calculation of the fragmentation term in the transversely polarized cross section
def get_dsigST(x, z, xF, pT, rs):

    Mh = {}
    Mh['pi+'] = 0.13957018  # self.aux.Mpi
    Mh['pi-'] = 0.13957018  # self.aux.Mpi
    Mh['pi0'] = 0.1349766  # self.aux.Mpi
    Mh['k+'] = 0.13957018  # self.aux.Mk
    Mh['k-'] = 0.13957018  # self.aux.Mk

    Mh = Mh[had]

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

    # Get arrays of the nonperturbative functions
    ft = get_ft(xp, Q2)
    h = get_h(x, Q2)
    Hxxpz = get_Hxxpz(z, Q2, s, t, u)

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

    ffcs += ftu * (hd * Hxxpz5d + hs * Hxxpz5s) + ftd * \
        (hu * Hxxpz5u + hs * Hxxpz5s) + fts * (hu * Hxxpz5u + hd * Hxxpz5d)

    ffcs += ftub * (hd * Hxxpz6d + hs * Hxxpz6s) + ftdb * \
        (hu * Hxxpz6u + hs * Hxxpz6s) + ftsb * (hu * Hxxpz6u + hd * Hxxpz6d)

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


def get_sig(xF, pT, rs, mode='gauss', nx=100, nz=100):
    if mode == 'gauss':
        dsigdzdx = np.vectorize(lambda x, z: get_dsig(x, z, xF, pT, rs))
        dsigdz = np.vectorize(lambda z: fixed_quad(
            lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
        sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
    elif mode == 'quad':
        sig = dblquad(lambda x, z: get_dsig(x, z, xF, pT, rs),
                      zmin, 1., xmin, lambda x: 1.)[0]
    return sig


def get_sigST(xF, pT, rs, mode='gauss', nx=100, nz=100):
    if mode == 'gauss':
        dsigdzdx = np.vectorize(lambda x, z: get_dsigST(x, z, xF, pT, rs))
        dsigdz = np.vectorize(lambda z: fixed_quad(
            lambda x: dsigdzdx(x, z), xmin(z), 1, n=nx)[0])
        sig = fixed_quad(dsigdz, zmin, 1, n=nz)[0]
    elif mode == 'quad':
        sig = dblquad(lambda x, z: get_dsigST(x, z, xF, pT, rs),
                      zmin, 1., xmin, lambda x: 1.)[0]
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

    # def test():
    #den = get_sig(xF,pT,rs,tar,had,mode='gauss',nx=100,nz=100)
    #num = get_sigST(xF,pT,rs,tar,had,mode='gauss',nx=100,nz=100)

    #AN = num/den
    #print AN

    # test()

#  def test2():
#    den = get_dsig(0.3,0.6,xF,pT,rs,tar,had)
#    num = get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
#
#    print den,num
#
#  from timeit import Timer
#  t = Timer("test2()", "from __main__ import test2")
#  print 't elapsed ',t.timeit(number=1)

    start = time.time()
    print get_dsig(0.3, 0.6, xF, pT, rs)
    print get_dsigST(0.3, 0.6, xF, pT, rs)
    end = time.time()
    print 'time=', (end - start)

    start = time.time()
    print get_dsig(0.3, 0.6, xF, pT, rs)
    print get_dsigST(0.3, 0.6, xF, pT, rs)
    end = time.time()
    print 'time=', (end - start)

    #print get_dsig(0.3,0.6,xF,pT,rs,tar,had)
    #print get_dsigST(0.3,0.6,xF,pT,rs,tar,had)

    # Integration of the numerator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)

    # Integration of the denominator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)
    #den = dblquad(lambda x,z: ANTHEORY().get_dsig(x,z,xF,pT,rs,tar,had),zmin,1.,xmin,lambda x: 1.)

    #AN = num[0]/den[0]
    # print(AN)
