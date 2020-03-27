import sys
import os
import numpy as np
from mpmath import fp
from scipy.special import gamma
from scipy.integrate import quad


class AUX:

    def __init__(self):

        self.set_constants()
        self.set_masses()
        self.set_couplings()
        self.set_evolution()
        self.set_ckm()

    def set_constants(self):

        self.CA = 3.0
        self.CF = 4.0 / 3.0
        self.TR = 0.5
        self.TF = 0.5
        self.euler = fp.euler

        # soft factors
        self.soft = {}

    def set_masses(self):

        self.me = 0.000511
        self.mmu = 0.105658
        self.mtau = 1.77684
        self.mu = 0.055
        self.md = 0.055
        self.ms = 0.2
        self.mc = 1.51
        self.mb = 4.92
        self.mZ = 91.1876
        self.mW = 80.398
        self.M = 0.93891897
        self.Mpi = 0.135
        self.Mk = 0.497
        #self.MA = 0.93891897 # Used for DY
        #self.MB = 0.93891897 # Used for DY

        self.me2 = self.me**2
        self.mmu2 = self.mmu**2
        self.mtau2 = self.mtau**2
        self.mu2 = self.mu**2
        self.md2 = self.md**2
        self.ms2 = self.ms**2
        self.mc2 = self.mc**2
        self.mb2 = self.mb**2
        self.mZ2 = self.mZ**2
        self.mW2 = self.mW**2
        self.M2 = self.M**2
        self.Mpi2 = self.Mpi**2
        
    def set_ckm(self): # CKM matrix
        self.Vud = 0.97420 # CKM matrix element +- 0.00021
        self.Vus = 0.2243  # CKM matrix element +- 0.0005
        self.Vcd = 0.218   # CKM matrix element +- 0.004

    def set_couplings(self):

        self.c2w = self.mW2 / self.mZ2
        self.s2w = 1.0 - self.c2w
        self.s2wMZ = 0.23116
        self.alfa = 1 / 137.036
        self.alphaSMZ = 0.118
        # axial couplings of quarks to Z boson
        self.cauz = 0.5 
        self.cacz = 0.5
        self.cadz = -0.5
        self.casz = -0.5
        # vector couplings of quarks to Z boson
        self.cvuz = self.cauz - 4./3.*self.s2w 
        self.cvcz = self.cacz - 4./3.*self.s2w
        self.cvdz = self.cadz + 2./3.*self.s2w
        self.cvsz = self.casz + 2./3.*self.s2w

    def set_evolution(self):
        self.Q02=2.0
        self.lam2=0.04

    def p2n(self, p):
        return np.copy(p[[0,3,4,1,2,5,6,7,8,9,10]])

    def charge_conj(self,p):
        return np.copy(p[[0,2,1,4,3,6,5,8,7,10,9]])

    def q2qbar(self, p):
        """
        distributions p
        p[0]    # g
        p[1]    # u
        p[2]    # ub
        p[3]    # d
        p[4]    # db
        p[5]    # s
        p[6]    # sb
        p[7]    # c
        p[8]    # cb
        p[9]    # b
        p[10]   # bb
        q <-> qbar used in Drell-Yan 
        """
        return np.copy(p[[0,2,1,4,3,6,5,8,7,10,9]])

    def _pip2piz(self, pip):
        """
        deprecated ??
        """
        piz = np.copy(pip)
        u = pip[1]
        ub = pip[2]
        d = pip[3]
        db = pip[4]
        piz[1] = 0.5 * (u + ub)
        piz[2] = 0.5 * (u + ub)
        piz[3] = 0.5 * (d + db)
        piz[4] = 0.5 * (d + db)
        return piz
  
    def _get_psi(self,i,N):
        return fp.psi(i,complex(N.real,N.imag))

    def get_psi(self,i,N):
        return np.array([self._get_psi(i,n) for n in N],dtype=complex)

    def beta(self,a,b):
        return gamma(a)*gamma(b)/gamma(a+b)


