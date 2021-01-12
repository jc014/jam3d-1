import sys
import os
import numpy as np
import time
from scipy.special import gamma
from tools.config import conf

class CORE:

    def beta(self, a, b):
        return gamma(a) * gamma(b) / gamma(a + b)

    def get_s(self,Q2):
        lam2 = conf['aux'].lam2
        Q02  = conf['aux'].Q02
        return np.log(np.log(Q2/lam2)/np.log(Q02/lam2))

    # shape funcs

    def __get_shape(self,x,p):
        return p[0]*x**p[1]*(1-x)**p[2]*(1+p[3]*x+p[4]*x**2)

    def _get_shape(self,x,p,s):
        N=p[0] + p[5] * s
        a=p[1] + p[6] * s
        b=p[2] + p[7] * s
        c=p[3] + p[8] * s
        d=p[4] + p[9] * s
        n= 1
        return self.__get_shape(x,[N/n,a,b,c,d])

    def get_shape(self,x,Q2,p1,p2):
        s=self.get_s(Q2)
        shape=self._get_shape(x,p1,s)
        shape2=self._get_shape(x,p2,s)
        if p2[0]!=0: shape=shape*(1+shape2)
        n1 = self.beta(p1[1]+2,p1[2]+1) + p1[3]*self.beta(p1[1]+3,p1[2]+1) + p1[4]*self.beta(p1[1]+4,p1[2]+1)
        n2 = p2[0] * (self.beta(p1[1]+p2[1]+2,p1[2]+p2[2]+1) + p1[4]*p2[4]*self.beta(p1[1]+p2[1]+6,p1[2]+p2[2]+1))
        n3 = p2[0] * (p1[3]+p2[3])*self.beta(p1[1]+p2[1]+3,p1[2]+p2[2]+1)
        n4 = p2[0] * (p1[3]*p2[4]+p2[3]*p1[4])*self.beta(p1[1]+p2[1]+5,p1[2]+p2[2]+1)
        n5 = p2[0] * (p1[3]*p2[3]+p1[4]+p2[4])*self.beta(p1[1]+p2[1]+4,p1[2]+p2[2]+1)
        n=n1+n2+n3+n4+n5
        shape=shape/n
        return shape

    def get_collinear(self,x,Q2):
        N = np.zeros(11)
        for i in range(11):
            N[i] = self.get_shape(x,Q2,self.shape1[i],self.shape2[i])
        return N

    # shape differentials

    def __get_dshape(self,x,p):
        return p[0]*x**p[1]*(1-x)**p[2]*(p[3]+2*p[4]*x+(1+p[3]*x+p[4]*x**2)*(p[1]/x-p[2]/(1-x)))

    def _get_dshape(self,x,p,s):
        N=p[0] + p[5] * s
        a=p[1] + p[6] * s
        b=p[2] + p[7] * s
        c=p[3] + p[8] * s
        d=p[4] + p[9] * s
        n=1
        return self.__get_dshape(x,[N/n,a,b,c,d])

    def get_dshape(self,x,Q2,p1,p2):
        s=self.get_s(Q2)
        dshape=self._get_dshape(x,p1,s)
        shape=self._get_shape(x,p1,s)
        dshape2=self._get_dshape(x,p2,s)
        shape2=self._get_shape(x,p2,s)
        if p2[0]!=0: dshape = shape*dshape2 + (1+shape2)*dshape
        n1 = self.beta(p1[1]+2,p1[2]+1) + p1[3]*self.beta(p1[1]+3,p1[2]+1) + p1[4]*self.beta(p1[1]+4,p1[2]+1)
        n2 = p2[0] * (self.beta(p1[1]+p2[1]+2,p1[2]+p2[2]+1) + p1[4]*p2[4]*self.beta(p1[1]+p2[1]+6,p1[2]+p2[2]+1))
        n3 = p2[0] * (p1[3]+p2[3])*self.beta(p1[1]+p2[1]+3,p1[2]+p2[2]+1)
        n4 = p2[0] * (p1[3]*p2[4]+p2[3]*p1[4])*self.beta(p1[1]+p2[1]+5,p1[2]+p2[2]+1)
        n5 = p2[0] * (p1[3]*p2[3]+p1[4]+p2[4])*self.beta(p1[1]+p2[1]+4,p1[2]+p2[2]+1)
        n=n1+n2+n3+n4+n5
        dshape=dshape/n
        return dshape

    def get_dcollinear(self,x,Q2):
        N = np.zeros(11)
        for i in range(11):
            N[i] = self.get_dshape(x,Q2,self.shape1[i],self.shape2[i])
        return N

    # moments

    def __get_mom(self,p,s,i):
        N=p[0] + p[5] * s
        a=p[1] + p[6] * s
        b=p[2] + p[7] * s
        c=p[3] + p[8] * s
        d=p[4] + p[9] * s
        ni= self.beta(i+p[1],p[2]+1) + p[3]*self.beta(i+p[1]+1,p[2] + 1) + p[4]*self.beta(i+p[1]+2,p[2] + 1)
        n1= self.beta(1+p[1],p[2]+1) + p[3]*self.beta(1+p[1]+1,p[2] + 1) + p[4]*self.beta(1+p[1]+2,p[2] + 1)
        return N*ni/n1

    def _get_mom(self,Q2,p1,p2,i):
        s=self.get_s(Q2)
        mom=self.__get_mom(p1,s,i)
        if p2[0]!=0: mom+=self.__get_mom(p2,s,i)
        return mom

    def get_mom(self,Q2,iflav,i):
        return self._get_mom(Q2,self.shape1[iflav],self.shape2[iflav],i)

    # widths

    def get_widths(self,Q2):
        s=np.log(Q2/conf['aux'].Q02 )
        return np.abs(self.widths1+s*self.widths2)

    # tmd

    def get_tmd(self,x,Q2,kT,lepton,dist,icol=False,deriv=False):

        col=self.get_C(x,Q2)

        if icol: return col

        widths=self.get_widths(Q2)
        if dist=='ffmu' or dist=='collinsmu': gauss=np.exp(-x**2 * kT**2/widths)
        else: gauss=np.exp(-kT**2/widths)

        if   lepton=='e' : M=conf['aux'].me
        elif lepton=='mu': M=conf['aux'].mmu
        #elif hadron=='k' : M=conf['aux'].Mk

        if   dist=='ldf':           norm=1/np.pi/widths
        elif dist=='transversity':  norm=1/np.pi/widths
        elif dist=='sivers':        norm=(2*M**2/widths) * (1/np.pi/widths)

        elif dist=='ffmu' :         norm=1/np.pi/widths
        elif dist=='collinsmu':     norm=2*x**2*M**2/widths * (1/np.pi/widths)

        return norm*col*gauss
