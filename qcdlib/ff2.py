#!/usr/bin/env python
import sys
import os
import numpy as np
import time
from qcdlib.core import CORE
from qcdlib.interpolator import INTERPOLATOR
from tools.config import conf

from scipy.special import gamma
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.dglap import DGLAP
from qcdlib.kernels import KERNELS
from qcdlib.mellin import MELLIN
from scipy.integrate import quad


class FF(CORE): 

    def __init__(self,spl):
  
        self.Q20=conf['Q20']
        self.mc2=conf['aux'].mc2
        self.mb2=conf['aux'].mb2
          
        self.mellin=conf['mellin']
        self.kernel=KERNELS(self.mellin,spl)
        self.dglap=DGLAP(self.mellin,conf['alphaS'],self.kernel,'truncated','LO')
  
        self.set_params()
        self.setup()
        self.ford=['g','u','ub','d','db','s','sb','c','cb','b','bb']

    def set_params(self):

        #--f(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
        params={}
        params['g1']   = np.array([0,-0.5,3,0,0])
        params['u1']   = np.array([1,-0.5,3,0,0])
        params['d1']   = np.array([1,-0.5,3,0,0])
        params['s1']   = np.array([1,-0.5,3,0,0])
        params['c1']   = np.array([1,-0.5,3,0,0])
        params['b1']   = np.array([1,-0.5,3,0,0])
        params['ub1']  = np.array([1,-0.5,3,0,0])
        params['db1']  = np.array([1,-0.5,3,0,0])
        params['sb1']  = np.array([1,-0.5,3,0,0])
        params['cb1']  = np.array([1,-0.5,3,0,0])
        params['bb1']  = np.array([1,-0.5,3,0,0])
        self.params=params

        #--widthds
        self._widths1_uv  = 0.3
        self._widths1_dv  = 0.3
        self._widths1_sea = 0.3
        self._widths2_uv  = 0
        self._widths2_dv  = 0
        self._widths2_sea = 0

        # internal
        self.widths1 = np.ones(11)
        self.widths2 = np.ones(11)
  
    def set_sumrules(self):
        pass  

    def set_moms(self):
        moms={}
        moms['g']  =self.get_moments('g1')
        moms['u1'] =self.get_moments('u1')
        moms['d1'] =self.get_moments('d1')
        moms['s1'] =self.get_moments('s1')
        moms['c1'] =self.get_moments('c1')
        moms['b1'] =self.get_moments('b1')
        moms['ub1']=self.get_moments('ub1')
        moms['db1']=self.get_moments('db1')
        moms['sb1']=self.get_moments('sb1')
        moms['cb1']=self.get_moments('cb1')
        moms['bb1']=self.get_moments('bb1')

        moms['up']=moms['u1']+moms['ub1']
        moms['dp']=moms['d1']+moms['db1']
        moms['sp']=moms['s1']+moms['sb1']
        moms['cp']=moms['c1']+moms['cb1']
        moms['bp']=moms['b1']+moms['bb1']
        moms['um']=moms['u1']-moms['ub1']
        moms['dm']=moms['d1']-moms['db1']
        moms['sm']=moms['s1']-moms['sb1']
        moms['cm']=moms['c1']-moms['cb1']
        moms['bm']=moms['b1']-moms['bb1']
  
        self.moms0=moms
        self.get_BC(moms)

    def set_widths(self):
        for i in range(11):
            if   i == 1: self.widths1[i] = self._widths1_uv
            elif i == 3: self.widths1[i] = self._widths1_dv
            else:        self.widths1[i] = self._widths1_sea
        for i in range(11):
            if   i == 1: self.widths2[i] = self._widths2_uv
            elif i == 3: self.widths2[i] = self._widths2_dv
            else:        self.widths2[i] = self._widths2_sea
  
    def setup(self):
        self.set_sumrules()
        self.set_moms()
        self.set_widths()
        #--store moments of a given Q2 that has been already calculated
        self.storage={}
  
    def beta(self,a,b):
        return gamma(a)*gamma(b)/gamma(a+b)
  
    def get_moments(self,flav,N=None):
        """
        if N==None: then parametrization is to be use to compute moments along mellin contour
        else the Nth moment is returned
        """
        if N==None: N=self.mellin.N
        M,a,b,c,d=self.params[flav]
        norm=self.beta(2+a,b+1)+c*self.beta(2+a+0.5,b+1)+d*self.beta(2+a+1.0,b+1)
        mom=self.beta(N+a,b+1)+c*self.beta(N+a+0.5,b+1)+d*self.beta(N+a+1.0,b+1)
        return M*mom/norm 
  
    def _get_BC(self,g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm):
        N=self.mellin.N
        
        # flav composition
        vm,vp={},{}
        vm[35]= bm + cm + dm + sm - 5*tm + um
        vm[24]= -4*bm + cm + dm + sm + um 
        vm[15]= -3*cm + dm + sm + um
        vm[8] = dm - 2*sp + 2*(-sm + sp) + um
        vm[3] = -dm + um
        vm[0] = np.zeros(N.size,dtype=complex)
        vp[0] = np.zeros(N.size,dtype=complex)
        vp[3] = -dp + up
        vp[8] = dp - 2*sp + up 
        vp[15]= -3*cp + dp + sp + up 
        vp[24]= -4*bp + cp + dp + sp + up 
        vp[35]= bp + cp + dp + sp - 5*tp + up 
        qs    = bp + cp + dp + sp + tp + up 
        qv    = bm + cm + dm + sm + tm + um 
        q     = np.zeros((2,N.size),dtype=complex)
        q[0]=np.copy(qs)
        q[1]=np.copy(g)
  
        BC={}
        BC['vm']=vm 
        BC['vp']=vp 
        BC['qv']=qv
        BC['q'] =q
        return BC
  
    def get_state(self):
        return (self.widths1,self.widths2,self.BC3,self.BC4,self.BC5)
    
    def set_state(self,state):
        self.widths1,self.widths2,self.BC3, self.BC4, self.BC5 = state[:]
        self.storage = {}

    def get_BC(self,moms):

        N=self.mellin.N
        zero=np.zeros(N.size,dtype=complex)

        ######################################
        # BC3 
        g   = moms['g']
        up  = moms['up']
        dp  = moms['dp']
        sp  = moms['sp']
        cp  = zero 
        bp  = zero
        um  = moms['um']
        dm  = moms['dm']
        sm  = moms['sm']
        cm  = zero 
        bm  = zero
        self.BC3=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,zero,zero)

        ######################################
        # BC for Nf=4
        BC4=self.dglap.evolve(self.BC3,self.Q20,self.mc2,3)
        g =BC4['g']
        up=BC4['up']
        dp=BC4['dp']
        sp=BC4['sp']
        cp=moms['cp']
        bp=zero
        um=BC4['um']
        dm=BC4['dm']
        sm=BC4['sm']
        cm=moms['cm']
        bm=zero
        self.BC4=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,zero,zero)

        ######################################
        # BC for Nf=5
        BC5=self.dglap.evolve(self.BC4,self.mc2,self.mb2,4)
        g =BC5['g']
        up=BC5['up']
        dp=BC5['dp']
        sp=BC5['sp']
        cp=BC5['cp']
        bp=moms['bp']
        um=BC5['um']
        dm=BC5['dm']
        sm=BC5['sm']
        cm=BC5['cm']
        bm=moms['bm']
        self.BC5=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,zero,zero)
  
    def evolve(self,Q2):
  
        if Q2 not in self.storage:
            if self.mb2<Q2: 
                self.storage[Q2]=self.dglap.evolve(self.BC5,self.mb2,Q2,5)
            elif self.mc2<=Q2 and Q2<=self.mb2: 
                self.storage[Q2]=self.dglap.evolve(self.BC4,self.mc2,Q2,4)
            elif Q2<self.mc2: 
                self.storage[Q2]=self.dglap.evolve(self.BC3,self.Q20,Q2,3)
    
    def get_xF(self,x,Q2,flav,evolve=True):
        if evolve: self.evolve(Q2)
        return x*self.mellin.invert(x,self.storage[Q2][flav])
  
    def get_xF0(self,x,flav):
        if   flav=='um': mom=self.moms0['um']
        elif flav=='dm': mom=self.moms0['dm']
        elif flav=='sm': mom=self.moms0['sm']
        return x*conf['mellin'].invert(x,mom)

    def get_C(self,x, Q2):
        self.evolve(Q2)
        return np.array([self.mellin.invert(x,self.storage[Q2][_]) for _ in self.ford])

    def get_dC(self, x, Q2):
        """
        print we need to think how to do this here
        """
        return 0#self.get_dcollinear(x,Q2)


if __name__ == '__main__':

    from qcdlib.aux import AUX
    from scipy.integrate import quad
    conf['order']='NLO'
    conf['Q20'] = 1.27**2
    conf['dglap mode']='truncated'
    conf['aux']=AUX()
    conf['mellin']=MELLIN(npts=16)
    conf['alphaS']=ALPHAS()
    conf['d1']  = FF('d1')
    

    x = 0.15
    Q2 = 200.4
    print conf['d1'].get_C(x, Q2)


      


