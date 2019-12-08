#!/usr/bin/env python
import pandas as pd
from fitlab.parman import PARMAN
from fitlab.mcsamp import MCSAMP
from fitlab.maxlike import ML
from tools.config import load_config, conf
from tools.tools import checkdir
from fitlab.resman import RESMAN

############################################################################
# params
conf['params']={}
# Parameters in gaussian approximation, parton model:
# TMD PDF:
conf['params']['pdf']={}
conf['params']['pdf']['widths1_uv']  ={'value':    5.77445e-01,'min': 0,'max':1,'fixed':False}
conf['params']['pdf']['widths2_uv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['pdf']['widths1_dv']  ={'value':    5.77445e-01,'min': 0,'max':1,'fixed':'widths1_uv'}
conf['params']['pdf']['widths2_dv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}
conf['params']['pdf']['widths1_sea'] ={'value':    5.87102e-01,'min': 0,'max':1,'fixed':False}
conf['params']['pdf']['widths2_sea'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}

# TMD FF:
conf['params']['ffpi']={}
conf['params']['ffpi']['widths1_fav']   ={'value':    1.16538e-01,'min': 0,'max':1,'fixed':False}
conf['params']['ffpi']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffpi']['widths1_ufav'] ={'value':    1.37436e-01,'min': 0,'max':1,'fixed':False}
conf['params']['ffpi']['widths2_ufav'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']={}
conf['params']['ffk']['widths1_fav']    ={'value':    1.32597e-01,'min': 0,'max':1,'fixed':False}
conf['params']['ffk']['widths2_fav']    ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']['widths1_ufav']  ={'value':    1.86177e-01,'min': 0,'max':1,'fixed':False}
conf['params']['ffk']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']={}
conf['params']['ffh']['widths1_fav']    ={'value':    1.32925e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_fav']    ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']['widths1_ufav']  ={'value':    1.86073e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}

############################################################################
# set data sets

conf['datasets']={}
conf['datasets']['sidis']={}

conf['datasets']['sidis']['xlsx']={}
conf['datasets']['sidis']['xlsx'][1000]='sidis/expdata/1000.xlsx'  # |  proton   | pi+   | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1001]='sidis/expdata/1001.xlsx'  # |  proton   | pi-   | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1004]='sidis/expdata/1004.xlsx'  # |  deuteron | pi+   | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx'  # |  deuteron | pi-   | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1002]='sidis/expdata/1002.xlsx'  # |  proton   | k+    | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1003]='sidis/expdata/1003.xlsx'  # |  proton   | k-    | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx'  # |  deuteron | k+    | M_Hermes | hermes 
conf['datasets']['sidis']['xlsx'][1007]='sidis/expdata/1007.xlsx'  # |  deuteron | k-    | M_Hermes | hermes 


############################################################################
# define datasets (kinematics) to perform the simulation

conf['datasets']={}
conf['datasets']['sidis']={}
conf['datasets']['sidis']['xlsx']={}
conf['datasets']['sidis']['xlsx'][1000]='./template.xlsx'  # |  proton   | pi+   
conf['datasets']['sidis']['norm']={}
for idx in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][idx]={'value':1,'fixed':True,'min':0,'max':1} 
conf['datasets']['sidis']['filters']={}


############################################################################
# run simulation
resman = RESMAN()
parman = conf['parman']
resman.get_residuals(parman.par)
print 


############################################################################
# prepare the new pseudo data set and save it as xlsx file
simdata=pd.DataFrame(conf['sidis tabs'][1000])
simdata['value']=simdata['thy']
simdata['stat_u']=pd.Series(simdata['thy']*simdata['%stat']/100,index=simdata.index)
simdata['syst_u']=pd.Series(simdata['thy']*simdata['%syst']/100,index=simdata.index)
simdata=simdata.drop(['r-residuals','shift','alpha','residuals','yp','yh','W2','Shift','N','thy','%stat','%syst'],axis=1)
print simdata[:10]

writer = pd.ExcelWriter('simulation.xlsx')
simdata.to_excel(writer,'Sheet1',index=False)
writer.save()



