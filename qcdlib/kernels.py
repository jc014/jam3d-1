import numpy as np
from mpmath import fp

class KERNELS:

    def __init__(self,mell,Type='f1'):

        self.D={}
        self.D['N']=mell.N
        self.D['nflav']=6+1
        self.D['norder']=2+1
        self.D['Nsize']=mell.N.size

        self.set_abbreviations()
        self.LO_unpolarized_splitting_functions()
        self.LO_polarized_splitting_functions()
        self.LO_unpol_timelike_splitting_functions()

        if   Type=='f1' : self.load_f1_spl()
        elif Type=='g1' : self.load_g1_spl()
        elif Type=='d1' : self.load_d1_spl()

    def set_abbreviations(self):
        D=self.D
        D['CA']=3.0
        D['CF']=4.0/3.0
        D['TR']=0.5
        zeta2=fp.zeta(2)
        zeta3=fp.zeta(3)
        N   = D['N']
        psi=lambda i,_N: fp.psi(i,complex(_N.real,_N.imag))
        S1f = lambda _N: fp.euler + psi(0,_N+1)
        S2f = lambda _N: zeta2 - psi(1,_N+1)
        S1 = np.array([S1f(n) for n in N])
        D['S1']=S1

    def LO_unpolarized_splitting_functions(self):

        D=self.D
        D['P0QQ']=np.zeros((D['nflav'],D['Nsize']),dtype=complex)
        D['P0QG']=np.zeros((D['nflav'],D['Nsize']),dtype=complex)
        D['P0GQ']=np.zeros((D['nflav'],D['Nsize']),dtype=complex)
        D['P0GG']=np.zeros((D['nflav'],D['Nsize']),dtype=complex)  
        D['P0'] = np.zeros((D['nflav'],2,2,D['Nsize']),dtype=complex)

        N=D['N']
        S1=D['S1']
        for Nf in range(3,D['nflav']):

            D['P0QQ'][Nf]=4.0/3.0*(3.0+2.0/N/(N+1)-4.0*S1)
            D['P0QG'][Nf]=2.0*(N**2+N+2)/(N*(N+1)*(N+2))*Nf
            D['P0GQ'][Nf]=8.0/3.0*(N**2+N+2)/(N-1)/N/(N+1)
            D['P0GG'][Nf]=3.0*(11.0/3.0+4.0/N/(N-1)+4.0/(N+1)/(N+2)-4.0*S1)-2.0/3.0*Nf

            D['P0'][Nf,0,0] = D['P0QQ'][Nf] 
            D['P0'][Nf,0,1] = D['P0QG'][Nf]
            D['P0'][Nf,1,0] = D['P0GQ'][Nf]
            D['P0'][Nf,1,1] = D['P0GG'][Nf]

    def LO_polarized_splitting_functions(self):
        D=self.D

        CF=D['CF']
        CA=D['CA']
        TR=D['TR']

        D['PP0'] = np.zeros((D['nflav'],2,2,D['Nsize']),dtype=complex)

        N=D['N']
        S1=D['S1']
        NM  = N - 1.
        N1  = N + 1
        N2  = N + 2
        PPQQA0 = 3. - 4.* S1 + 2./(N * N1)
        PPQGA0 = 4.* NM / (N * N1)
        PPGQA0 = 2.* N2 / (N * N1)
        PPGGA0 = 11./3. - 4.* S1 + 8./ (N * N1)
        PPGGB0 = - 4./3.

        for NF in range(3,D['nflav']):

            D['PP0'][NF,0,0]=CF*PPQQA0
            D['PP0'][NF,0,1]=TR*NF*PPQGA0
            D['PP0'][NF,1,0]=CF*PPGQA0
            D['PP0'][NF,1,1]=CA*PPGGA0+TR*NF*PPGGB0
          
    def LO_unpol_timelike_splitting_functions(self):
      
        D=self.D
        D['PT0QQ']=np.zeros((D['nflav'],D['Nsize']),dtype=complex) 
        D['PT0GQ']=np.zeros((D['nflav'],D['Nsize']),dtype=complex) 
        D['PT0QG']=np.zeros((D['nflav'],D['Nsize']),dtype=complex) 
        D['PT0GG']=np.zeros((D['nflav'],D['Nsize']),dtype=complex) 
        D['PT0'] = np.zeros((D['nflav'],2,2,D['Nsize']),dtype=complex)

        N=D['N']
        S1=D['S1']
        for Nf in range(3,D['nflav']):

            D['P0QQ'][Nf]=4.0/3.0*(3.0+2.0/N/(N+1)-4.0*S1)
            D['P0QG'][Nf]=2.0*(N**2+N+2)/(N*(N+1)*(N+2))*Nf
            D['P0GQ'][Nf]=8.0/3.0*(N**2+N+2)/(N-1)/N/(N+1)
            D['P0GG'][Nf]=3.0*(11.0/3.0+4.0/N/(N-1)+4.0/(N+1)/(N+2)-4.0*S1)-2.0/3.0*Nf

            D['PT0QQ'][Nf]=D['P0QQ'][Nf]
            D['PT0QG'][Nf]=D['P0QG'][Nf]/2./Nf
            D['PT0GQ'][Nf]=D['P0GQ'][Nf]*2.*Nf
            D['PT0GG'][Nf]=D['P0GG'][Nf]
        
            D['PT0'][Nf,0,0] = D['PT0QQ'][Nf]
            D['PT0'][Nf,0,1] = D['PT0GQ'][Nf]
            D['PT0'][Nf,1,0] = D['PT0QG'][Nf]
            D['PT0'][Nf,1,1] = D['PT0GG'][Nf]

    def load_f1_spl(self):
        D=self.D
        Nsize=D['N'].size
        nflav=D['nflav']
        norder=D['norder']
 
        # initialize flav composed splitting functions arrays
        self.PNSP=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.PNSM=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.PNSV=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.P   =np.zeros((nflav,norder,2,2,Nsize),dtype=complex)

        for Nf in range(3,nflav):
            self.PNSP[Nf,0] = D['P0QQ'][Nf] 
            self.PNSM[Nf,0] = D['P0QQ'][Nf] 
            self.PNSV[Nf,0] = D['P0QQ'][Nf] 
            self.P[Nf,0]    = D['P0'][Nf]

    def load_g1_spl(self):
        D=self.D
        Nsize=D['N'].size
        nflav=D['nflav']
        norder=D['norder']
 
        # initialize flav composed splitting functions arrays
        self.PNSP=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.PNSM=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.PNSV=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.P    =np.zeros((nflav,norder,2,2,Nsize),dtype=complex)

        for Nf in range(3,nflav):
            self.PNSP[Nf,0] = D['P0QQ'][Nf]
            self.PNSM[Nf,0] = D['P0QQ'][Nf]
            self.PNSV[Nf,0] = D['P0QQ'][Nf]
            self.P[Nf,0]    = D['PP0'][Nf]
          
    def load_d1_spl(self):
      
        D=self.D
        Nsize=D['N'].size
        nflav=D['nflav']
        norder=D['norder']
 
        # initialize flav composed splitting functions arrays
        self.PNSP=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.PNSM=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.PNSV=np.zeros((nflav,norder,Nsize),dtype=complex)
        self.P   =np.zeros((nflav,norder,2,2,Nsize),dtype=complex)

        for Nf in range(3,nflav):
 
            # LO unpolarized
            self.PNSP[Nf,0] = D['PT0QQ'][Nf] 
            self.PNSM[Nf,0] = D['PT0QQ'][Nf] 
            self.PNSV[Nf,0] = D['PT0QQ'][Nf] 
            self.P[Nf,0]    = D['PT0'][Nf]


