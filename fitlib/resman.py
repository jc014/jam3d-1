#!/usr/bin/env python
import sys,os
import time
import numpy as np

#--from qcdlib
import qcdlib
from   qcdlib import pdf0,ff0,pdf1,ff1,pdf2,ff2
import qcdlib.aux
import qcdlib.alphaS
import qcdlib.interpolator
import qcdlib.alphaS
import qcdlib.mellin

#--from obslib
import obslib.sidis.residuals
import obslib.sidis.reader
import obslib.sia.residuals
import obslib.sia.reader
import obslib.moments.reader
import obslib.moments.residuals
import obslib.AN_pp.residuals
import obslib.AN_pp.reader
import obslib.dy.reader
import obslib.dy.residuals
import obslib.wz.reader
import obslib.wz.residuals

#--from fitlib
from fitlib.parman import PARMAN

#--from tools 
from tools.tools    import checkdir
from tools.config   import conf,load_config
from tools.parallel import PARALLEL

class RESMAN:

    def __init__(self,nworkers=2,parallel=True,datasets=True):

        self.setup_core()
        self.parman=PARMAN()
        if datasets:

            if 'sidis'   in conf['datasets']: self.setup_sidis()
            if 'sia'     in conf['datasets']: self.setup_sia()
            if 'moments' in conf['datasets']: self.setup_moments()
            if 'AN'      in conf['datasets']: self.setup_AN()
            if 'dy'      in conf['datasets']: self.setup_dy()
            if 'wz'      in conf['datasets']: self.setup_wz()

        if  parallel:
            self.setup_parallel(nworkers)
            self.requests=self.get_requests()
   
    def setup_core(self):

        conf['aux'] = qcdlib.aux.AUX()
        conf['mellin']= qcdlib.mellin.MELLIN(npts=4)
        conf['alphaS']= qcdlib.alphaS.ALPHAS()

        if 'pdf parametrization' in conf:

            if conf['pdf parametrization']==0: conf['pdf']= pdf0.PDF()
            if conf['pdf parametrization']==1: conf['pdf']= pdf1.PDF()
            if conf['pdf parametrization']==2: conf['pdf']= pdf2.PDF()
            if conf['pdf parametrization']==3: conf['pdf']= pdf3.PDF()

        if 'pdf'           in conf['params']: conf['pdf']          = pdf0.PDF()
        #if 'pdfpi-'        in conf['params']: conf['pdfpi-']       = pdf0.PDF('pi-')
        if 'transversity'  in conf['params']: conf['transversity'] = pdf1.PDF()
        if 'sivers'        in conf['params']: conf['sivers']       = pdf1.PDF()
        if 'boermulders'   in conf['params']: conf['boermulders']  = pdf1.PDF()
        if 'ffpi'          in conf['params']: conf['ffpi']         = ff0.FF('pi')
        if 'ffk'           in conf['params']: conf['ffk']          = ff0.FF('k')
        if 'collinspi'     in conf['params']: conf['collinspi']    = ff1.FF('pi')
        if 'collinsk'      in conf['params']: conf['collinsk']     = ff1.FF('k')
        if 'Htildepi'      in conf['params']: conf['Htildepi']     = ff1.FF('pi')
        if 'Htildek'       in conf['params']: conf['Htildek']      = ff1.FF('k')

        if 'transversity+' in conf['params']: conf['transversity'] = pdf2.PDF('h1') # Transversity  
        if 'collinspi+'    in conf['params']: conf['collinspi']    = ff2.FF('Col')  # Collins
        if 'sivers+'       in conf['params']: conf['sivers']    = pdf2.PDF('Siv')    # Sivers

    def setup_sidis(self):
        conf['sidis tabs']    = obslib.sidis.reader.READER().load_data_sets('sidis')
        self.sidisres = obslib.sidis.residuals.RESIDUALS()

    def setup_sia(self):
        conf['sia tabs']    = obslib.sia.reader.READER().load_data_sets('sia')
        self.siares = obslib.sia.residuals.RESIDUALS()

    def setup_AN(self):
        conf['AN tabs']   = obslib.AN_pp.reader.READER().load_data_sets('AN')
        self.ANres = obslib.AN_pp.residuals.RESIDUALS()
    
    def setup_dy(self):
        conf['dy tabs']   = obslib.dy.reader.READER().load_data_sets('dy')
        self.dyres = obslib.dy.residuals.RESIDUALS()

    def setup_wz(self):
        conf['wz tabs']   = obslib.wz.reader.READER().load_data_sets('wz')
        self.wzres = obslib.wz.residuals.RESIDUALS()

    def setup_parallel(self,nworkers):
        self.parallel=PARALLEL()
        self.parallel.task=self.task
        self.parallel.set_state=self.set_state
        self.parallel.setup_master()
        self.parallel.setup_workers(nworkers)
        self.nworkers=nworkers

    def get_state(self):
        state={}
        if 'pdf'          in conf: state['pdf'         ]    = conf['pdf'          ].get_state()
        if 'pdfpi-'       in conf: state['pdfpi-'      ]    = conf['pdfpi-'       ].get_state()
        if 'transversity' in conf: state['transversity']    = conf['transversity' ].get_state()
        if 'sivers'       in conf: state['sivers'      ]    = conf['sivers'       ].get_state()
        if 'boermulders'  in conf: state['boermulders' ]    = conf['boermulders'  ].get_state()
        if 'ffpi'         in conf: state['ffpi'        ]    = conf['ffpi'         ].get_state()
        if 'ffk'          in conf: state['ffk'         ]    = conf['ffk'          ].get_state()
        if 'collinspi'    in conf: state['collinspi'   ]    = conf['collinspi'    ].get_state()
        if 'collinsk'     in conf: state['collinsk'    ]    = conf['collinsk'     ].get_state()
        if 'Htildepi'     in conf: state['Htildepi'    ]    = conf['Htildepi'     ].get_state()
        if 'Htildek'      in conf: state['Htildek'     ]    = conf['Htildek'      ].get_state()
        return state

    def set_state(self,state):
        if 'pdf'          in conf: conf['pdf'         ].set_state(state['pdf'         ])
        if 'pdfpi-'       in conf: conf['pdfpi-'      ].set_state(state['pdfpi-'      ])
        if 'transversity' in conf: conf['transversity'].set_state(state['transversity'])
        if 'sivers'       in conf: conf['sivers'      ].set_state(state['sivers'      ])
        if 'boermulders'  in conf: conf['boermulders' ].set_state(state['boermulders' ])
        if 'ffpi'         in conf: conf['ffpi'        ].set_state(state['ffpi'        ])
        if 'ffk'          in conf: conf['ffk'         ].set_state(state['ffk'         ])
        if 'collinspi'    in conf: conf['collinspi'   ].set_state(state['collinspi'   ])
        if 'collinsk'     in conf: conf['collinsk'    ].set_state(state['collinsk'    ])
        if 'Htildepi'     in conf: conf['Htildepi'    ].set_state(state['Htildepi'    ])
        if 'Htildek'      in conf: conf['Htildek'     ].set_state(state['Htildek'     ])
  
    def distribute_requests(self,container,requests):
        cnt=0
        for request in requests:
            container[cnt].append(request)
            cnt+=1
            if cnt==self.nworkers: cnt=0

    def get_requests(self):
        container=[[] for _ in range(self.nworkers)]
        if 'sidis'  in conf['datasets']:  self.distribute_requests(container,self.sidisres.requests) 
        if 'sia'    in conf['datasets']:  self.distribute_requests(container,self.siares.requests) 
        if 'AN'     in conf['datasets']:  self.distribute_requests(container,self.ANres.requests) 
        if 'dy'     in conf['datasets']:  self.distribute_requests(container,self.dyres.requests)
        if 'wz'     in conf['datasets']:  self.distribute_requests(container,self.wzres.requests)
        return container

    def task(self,request):
        for i in range(len(request)):
            if  request[i]['reaction']=='sidis' :  self.sidisres.process_request(request[i])
            if  request[i]['reaction']=='sia'   :  self.siares.process_request(request[i])
            if  request[i]['reaction']=='AN'    :  self.ANres.process_request(request[i])
            if  request[i]['reaction']=='dy'    :  self.dyres.process_request(request[i])
            if  request[i]['reaction']=='wz'    :  self.wzres.process_request(request[i])
        return request
 
    def get_residuals(self,par):
        self.parman.set_new_params(par)
        state=self.get_state()
        self.parallel.update_workers(state)
        results=self.parallel.send_tasks(self.requests)

        #--update tables with the new theory values
        for chunk in results:
            for request in chunk:
                if request['reaction']=='sidis'  : self.sidisres.update_tabs_external(request)
                if request['reaction']=='sia'    : self.siares.update_tabs_external(request)
                if request['reaction']=='AN'     : self.ANres.update_tabs_external(request)
                if request['reaction']=='dy'     : self.dyres.update_tabs_external(request)
                if request['reaction']=='wz'     : self.wzres.update_tabs_external(request)

        #--compute residuals
        res,rres,nres=[],[],[]
        if 'sidis' in conf['datasets']:
            out=self.sidisres.get_residuals(calc=False)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'sia' in conf['datasets']:
            out=self.siares.get_residuals(calc=False)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'AN' in conf['datasets']:
            out=self.ANres.get_residuals(calc=False)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'dy' in conf['datasets']:
            out=self.dyres.get_residuals(calc=False)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'wz' in conf['datasets']:
            out=self.wzres.get_residuals(calc=False)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        return res,rres,nres

    def get_data_info(self):

        #--compute residuals
        reaction=[]
        if 'sidis' in conf['datasets']:
            out=self.sidisres.get_residuals(calc=False)
            reaction.extend(['sidis' for _ in out[0]])
        if 'sia' in conf['datasets']:
            out=self.siares.get_residuals(calc=False)
            reaction.extend(['sia' for _ in out[0]])
        if 'AN' in conf['datasets']:
            out=self.ANres.get_residuals(calc=False)
            reaction.extend(['AN' for _ in out[0]])
        if 'dy' in conf['datasets']:
            out=self.dyres.get_residuals(calc=False)
            reaction.extend(['dy' for _ in out[0]])
        if 'wz' in conf['datasets']:
            out=self.dyres.get_residuals(calc=False)
            reaction.extend(['wz' for _ in out[0]])
        return reaction

    def gen_report(self,verb=0,level=0):
        L=[]
        if 'sidis'   in conf['datasets']: L.extend(self.sidisres.gen_report(verb,level))
        if 'sia'     in conf['datasets']: L.extend(self.siares.gen_report(verb,level))
        if 'AN'      in conf['datasets']: L.extend(self.ANres.gen_report(verb,level))
        if 'dy'      in conf['datasets']: L.extend(self.dyres.gen_report(verb,level))
        if 'wz'      in conf['datasets']: L.extend(self.wzres.gen_report(verb,level))
        return L

    def get_chi2(self):
        data={}
        if 'sidis'   in conf['datasets']: data.update(self.sidisres.get_chi2())
        if 'sia'     in conf['datasets']: data.update(self.siares.get_chi2())
        if 'AN'      in conf['datasets']: data.update(self.ANres.get_chi2())
        if 'dy'      in conf['datasets']: data.update(self.dyres.get_chi2())
        if 'wz'      in conf['datasets']: data.update(self.wzres.get_chi2())
        return data

    def test(self,ntasks=10):
        #--loop over states 
        print '='*20
        t=time.time() 
        for _ in range(ntasks): 
            par=self.parman.par
            par*=(1+0.01*np.random.randn(par.size))
            res,rres,nres=self.get_residuals(par)
            chi2=np.sum(res**2)
            print '(%d/%d) chi2=%f'%(_,ntasks,chi2) 
        print '='*20
        elapsed_time=time.time()-t 
        print 'elapsed time :%f'%elapsed_time
        return elapsed_time

    def shutdown(self):
        self.parallel.stop_workers()

if __name__=='__main__':

    load_config('input_dglap.py')
    nworkers=3
    resman=RESMAN(nworkers)
    resman.test()
    resman.shutdown()




