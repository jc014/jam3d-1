conf={}

#--setups
conf['bootstrap']=True
conf['flat par']=True
conf['ftol']=1e-8
conf['ncpus']=4
conf['Q20']   = 1.27**2
conf['order'] = 'LO'

#--datasets

conf['datasets']={}

conf['datasets']['sidis']={}
conf['datasets']['sidis']['filters']=["z>0.2 and z<0.6 and Q2>1.63 and pT>0.2 and pT<0.9"]
conf['datasets']['sidis']['xlsx']={}

#--Unpol multiplicities (HERMES only)
conf['datasets']['sidis']['xlsx'][1000]='sidis/expdata/1000.xlsx'  # |  proton   | pi+   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1001]='sidis/expdata/1001.xlsx'  # |  proton   | pi-   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1002]='sidis/expdata/1002.xlsx'  # |  proton | k+   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1003]='sidis/expdata/1003.xlsx'  # |  proton | k-   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1004]='sidis/expdata/1004.xlsx'  # |  deuteron | pi+   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx'  # |  deuteron | pi-   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx'  # |  deuteron | k+   | M_Hermes | HERMES
conf['datasets']['sidis']['xlsx'][1007]='sidis/expdata/1007.xlsx'  # |  deuteron | k-   | M_Hermes | HERMES

#--SIDIS cos2phi
conf['datasets']['sidis']['xlsx'][5002]='sidis/expdata/5002.xlsx'  # |  proton   | pi-   | AUUcos2 | HERMES
conf['datasets']['sidis']['xlsx'][5003]='sidis/expdata/5003.xlsx'  # |  proton   | pi+   | AUUcos2 | HERMES
conf['datasets']['sidis']['xlsx'][5006]='sidis/expdata/5006.xlsx'  # |  deuteron | pi-   | AUUcos2 | HERMES
conf['datasets']['sidis']['xlsx'][5007]='sidis/expdata/5007.xlsx'  # |  deuteron | pi+   | AUUcos2 | HERMES


conf['datasets']['sidis']['norm']={}
for k in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][k]={'value':1.00000e+00,'fixed':False,'min':8.00000e-01,'max':1.20000e+00}


#--Collins effect SIA
conf['datasets']['sia']={}
conf['datasets']['sia']['filters']=[]
#conf['datasets']['sia']['filters'].append("pT/z1<3.5")
#conf['datasets']['sia']['filters'].append("z2<0.7")
conf['datasets']['sia']['xlsx']={}

conf['datasets']['sia']['xlsx'][1000]='sia/expdata/1000.xlsx' # babar | pi,pi | AUL-0     | 9      | z1,z2,pT0  |
conf['datasets']['sia']['xlsx'][1001]='sia/expdata/1001.xlsx' # babar | pi,pi | AUC-0     | 9      | z1,z2,pT0  |
conf['datasets']['sia']['xlsx'][1002]='sia/expdata/1002.xlsx' # babar | pi,pi | AUC-0     | 36     | z1,z2      |
conf['datasets']['sia']['xlsx'][1003]='sia/expdata/1003.xlsx' # babar | pi,pi | AUL-0     | 36     | z1,z2      |
conf['datasets']['sia']['xlsx'][1004]='sia/expdata/1004.xlsx' # belle | pi,pi | AUT-0-CCP | 16     | z1,z2,qT   |
conf['datasets']['sia']['xlsx'][1005]='sia/expdata/1005.xlsx' # belle | pi,pi | AUT-0     | 16     | z1,z2,qT   |
conf['datasets']['sia']['xlsx'][2008]='sia/expdata/2008.xlsx' # babar | pi,pi | AUL-0     | 16     | z1,z2      |
conf['datasets']['sia']['xlsx'][2009]='sia/expdata/2009.xlsx' # babar | pi,pi | AUC-0     | 16     | z1,z2      |
conf['datasets']['sia']['xlsx'][3000]='sia/expdata/3000.xlsx' # bes3 | pi,pi | AUL-0     | 6     | z1,z2        |
conf['datasets']['sia']['xlsx'][3001]='sia/expdata/3001.xlsx' # bes3 | pi,pi | AUC-0     | 6     | z1,z2        |
conf['datasets']['sia']['xlsx'][3002]='sia/expdata/3002.xlsx' # bes3 | pi,pi | AUL-0     | 5     | z1,z2,pT     |
conf['datasets']['sia']['xlsx'][3003]='sia/expdata/3003.xlsx' # bes3 | pi,pi | AUC-0     | 5     | z1,z2.pT     |

conf['datasets']['sia']['norm']={}
for k in conf['datasets']['sia']['xlsx']: conf['datasets']['sia']['norm'][k]={'value':1.00000e+00,'fixed':False,'min':8.00000e-01,'max':1.20000e+00}


#--parameters
conf['params']={}


#--Parameters in gaussian approximation, parton model:
#--TMD PDF:
conf['params']['pdf']={}
conf['params']['pdf']['widths1_uv']  ={'value':    5.51292e-01,'min': 0.1,'max':0.8,'fixed':False}
conf['params']['pdf']['widths2_uv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['pdf']['widths1_dv']  ={'value':    5.40344e-01,'min': 0.1,'max':0.8,'fixed':'widths1_uv'}
conf['params']['pdf']['widths2_dv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}
conf['params']['pdf']['widths1_sea'] ={'value':    5.43879e-01,'min': 0.1,'max':0.8,'fixed':False}
conf['params']['pdf']['widths2_sea'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}

#--TMD FF:
conf['params']['ffpi']={}
conf['params']['ffpi']['widths1_fav']  ={'value':    1.22011e-01,'min': 0.05,'max':0.3,'fixed':False}
conf['params']['ffpi']['widths2_fav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffpi']['widths1_ufav'] ={'value':    1.43665e-01,'min': 0.05,'max':0.3,'fixed':False}
conf['params']['ffpi']['widths2_ufav'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}

conf['params']['ffk']={}
conf['params']['ffk']['widths1_fav']   ={'value':    1.32333e-01,'min': 0.05,'max':0.3,'fixed':False}
conf['params']['ffk']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']['widths1_ufav']  ={'value':    2.02660e-01,'min': 0.05,'max':0.3,'fixed':False}
conf['params']['ffk']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}

conf['params']['ffh']={}
conf['params']['ffh']['widths1_fav']   ={'value':    1.24050e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']['widths1_ufav']  ={'value':    1.43730e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}

#Boer-Mulders function
conf['params']['boermulders']={}
conf['params']['boermulders']['widths1_uv']  ={'value':    2.47890e-01, 'min':0.1, 'max':1.0, 'fixed': False}
conf['params']['boermulders']['widths1_dv']  ={'value':    3.27437e-01, 'min':0, 'max':2, 'fixed': 'widths1_uv'}
conf['params']['boermulders']['widths1_sea'] ={'value':    4.88407e-01, 'min':0.1, 'max':2.0, 'fixed': True}
conf['params']['boermulders']['widths2_uv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': True}
conf['params']['boermulders']['widths2_dv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}
conf['params']['boermulders']['widths2_sea'] ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}

conf['params']['boermulders']['u N0 1'] ={'value':   -0.023648e+00, 'min': -0.15, 'max': 0, 'fixed': False}
conf['params']['boermulders']['u N1 1'] ={'value':    0.005e+00, 'min': 0.001, 'max': 0.005, 'fixed':False}
conf['params']['boermulders']['u a0 1'] ={'value':    7.58691e-01, 'min': -1.0, 'max':2.0, 'fixed': False}
conf['params']['boermulders']['u a1 1'] ={'value':    0.0000e+00, 'min':0, 'max':1.0, 'fixed':True}
conf['params']['boermulders']['u b0 1'] ={'value':    1.16529e+01, 'min': 0.5, 'max':17.0, 'fixed': False}
conf['params']['boermulders']['u b1 1'] ={'value':    1.0000e+00, 'min': 0.5, 'max':1.5, 'fixed':False}
conf['params']['boermulders']['u c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['u a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['u b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['u b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['u c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['u d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['boermulders']['d N0 1'] ={'value':    -0.01373e+00, 'min': -0.15, 'max': 0, 'fixed': False}
conf['params']['boermulders']['d N1 1'] ={'value':    -0.005e+00, 'min': -0.025, 'max': -0.005, 'fixed':False}
conf['params']['boermulders']['d a0 1'] ={'value':    1.19070e+00, 'min': -1.0, 'max':2, 'fixed': False}
conf['params']['boermulders']['d a1 1'] ={'value':    0.00000e+00, 'min':0, 'max':1.0, 'fixed': 'u a1 1'}
conf['params']['boermulders']['d b0 1'] ={'value':    1.95150e+01, 'min': 0, 'max':17.0, 'fixed': 'u b0 1'}
conf['params']['boermulders']['d b1 1'] ={'value':    0.50000e+00, 'min': 0.5, 'max':1.0, 'fixed': 'u b1 1'}
conf['params']['boermulders']['d c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['d a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['d b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['d b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['d c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['d d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['boermulders']['s N0 1'] ={'value':    0.00000e+00, 'min': -0.05, 'max': 0.05, 'fixed': True}
conf['params']['boermulders']['s N1 1'] ={'value':    0.00000e+00, 'min': -1.0, 'max': 0.0, 'fixed': True}
conf['params']['boermulders']['s a0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max':10, 'fixed': True}
conf['params']['boermulders']['s a1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max':0.0, 'fixed': True}
conf['params']['boermulders']['s b0 1'] ={'value':    15.00000e+00, 'min': 0, 'max':20.0, 'fixed': 'u b0 1'}
conf['params']['boermulders']['s b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':1.0, 'fixed': 'u b1 1'}
conf['params']['boermulders']['s c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['s a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['s b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['s b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['boermulders']['s c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['boermulders']['s d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['boermulders']['ub N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['boermulders']['ub N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['boermulders']['ub a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['boermulders']['ub a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['boermulders']['ub b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['boermulders']['ub b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['boermulders']['ub c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['boermulders']['ub c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['boermulders']['ub d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['boermulders']['ub d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['boermulders']['ub N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['boermulders']['ub N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['boermulders']['ub a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['boermulders']['ub a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['boermulders']['ub b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['boermulders']['ub b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['boermulders']['ub c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['boermulders']['ub c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['boermulders']['ub d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['boermulders']['ub d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

conf['params']['boermulders']['db N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['boermulders']['db N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['boermulders']['db a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['boermulders']['db a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['boermulders']['db b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['boermulders']['db b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['boermulders']['db c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['boermulders']['db c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['boermulders']['db d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['boermulders']['db d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['boermulders']['db N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['boermulders']['db N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['boermulders']['db a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['boermulders']['db a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['boermulders']['db b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['boermulders']['db b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['boermulders']['db c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['boermulders']['db c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['boermulders']['db d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['boermulders']['db d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

conf['params']['boermulders']['sb N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['boermulders']['sb N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['boermulders']['sb a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['boermulders']['sb a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['boermulders']['sb b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['boermulders']['sb b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['boermulders']['sb c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['boermulders']['sb c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['boermulders']['sb d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['boermulders']['sb d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['boermulders']['sb N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['boermulders']['sb N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['boermulders']['sb a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['boermulders']['sb a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['boermulders']['sb b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['boermulders']['sb b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['boermulders']['sb c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['boermulders']['sb c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['boermulders']['sb d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['boermulders']['sb d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}


# Collins TMD FF:
conf['params']['collinspi']={}

conf['params']['collinspi']['widths1_fav']   ={'value':    1.22288e-01,'min':0.025, 'max':0.25, 'fixed':False}
conf['params']['collinspi']['widths1_ufav'] ={'value':    6.11536e-02,'min':0.01, 'max':0.15, 'fixed':False}
conf['params']['collinspi']['widths2_fav']   ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['collinspi']['widths2_ufav'] ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':'widths2_fav'}

conf['params']['collinspi']['u N0 1'] ={'value':    0.2e+00,'min': -0.05, 'max':0.5, 'fixed':False}
conf['params']['collinspi']['u N1 1'] ={'value':    0e+00,'min': -0.0015, 'max': 0, 'fixed':False}
conf['params']['collinspi']['u a0 1'] ={'value':    2.5,'min':-1, 'max':7, 'fixed':False}
conf['params']['collinspi']['u a1 1'] ={'value':    0.00000e+00,'min':0, 'max': 0.25, 'fixed':False}
conf['params']['collinspi']['u b0 1'] ={'value':    2.5e+00,'min': 0.5, 'max': 5, 'fixed':False}
conf['params']['collinspi']['u b1 1'] ={'value':    0.0e+00,'min': 0.0, 'max':1, 'fixed':False}
conf['params']['collinspi']['u c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u N0 2'] ={'value':    1000,'min': 0, 'max': 1000, 'fixed':True}
conf['params']['collinspi']['u N1 2'] ={'value':    0.00000e+00,'min': -2, 'max': 0, 'fixed':True}
conf['params']['collinspi']['u a0 2'] ={'value':    0,'min':-0.1, 'max': 0.1, 'fixed':True}
conf['params']['collinspi']['u a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 0, 'fixed':True}
conf['params']['collinspi']['u b0 2'] ={'value':    11,'min':10, 'max':15, 'fixed':False}
conf['params']['collinspi']['u b1 2'] ={'value':    0.00e+00,'min': -0.5, 'max':0.5, 'fixed':False}
conf['params']['collinspi']['u c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['collinspi']['d N0 1'] ={'value':   -0.4,'min': -0.75, 'max':0.05, 'fixed':False}
conf['params']['collinspi']['d N1 1'] ={'value':    0e+00,'min': 0, 'max': 0.0015, 'fixed':False}
conf['params']['collinspi']['d a0 1'] ={'value':    0.5e+00,'min':0.5, 'max': 3, 'fixed':False}
conf['params']['collinspi']['d a1 1'] ={'value':    0.00000e+00,'min':0, 'max':0.25, 'fixed':False}
conf['params']['collinspi']['d b0 1'] ={'value':    3.41000e+00,'min': 2, 'max':5, 'fixed':False}
conf['params']['collinspi']['d b1 1'] ={'value':    0.000e+00,'min': 0, 'max':1, 'fixed':False}
conf['params']['collinspi']['d c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d N0 2'] ={'value':    200,'min':0, 'max':200, 'fixed':True}
conf['params']['collinspi']['d N1 2'] ={'value':    0.00000e+00,'min': -2, 'max': 0, 'fixed':True}
conf['params']['collinspi']['d a0 2'] ={'value':    0,'min': -0.1, 'max': 0.1, 'fixed':True}
conf['params']['collinspi']['d a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 0, 'fixed':True}
conf['params']['collinspi']['d b0 2'] ={'value':    10.80e+00,'min': 5, 'max':25, 'fixed':False}
conf['params']['collinspi']['d b1 2'] ={'value':    0.00000e+00,'min': -0.5, 'max':0.5, 'fixed':False}
conf['params']['collinspi']['d c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['collinspi']['s N0 1'] ={'value':   -5.94721e+00,'min': -15, 'max': 10, 'fixed':'d N0 1'}
conf['params']['collinspi']['s N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 1'}
conf['params']['collinspi']['s a0 1'] ={'value':    4.25846e+00,'min':-1, 'max': 20, 'fixed':'d a0 1'}
conf['params']['collinspi']['s a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 1'}
conf['params']['collinspi']['s b0 1'] ={'value':    4.19000e+00,'min': 2.19, 'max':10, 'fixed':'d b0 1'}
conf['params']['collinspi']['s b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b1 1'}
conf['params']['collinspi']['s c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c0 1'}
conf['params']['collinspi']['s c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c1 1'}
conf['params']['collinspi']['s d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d0 1'}
conf['params']['collinspi']['s d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d1 1'}
conf['params']['collinspi']['s N0 2'] ={'value':   -3.06065e+02,'min': 0, 'max': 1, 'fixed':'d N0 2'}
conf['params']['collinspi']['s N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 2'}
conf['params']['collinspi']['s a0 2'] ={'value':    8.69912e-01,'min':-2, 'max': 2, 'fixed':'d a0 2'}
conf['params']['collinspi']['s a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 2'}
conf['params']['collinspi']['s b0 2'] ={'value':    1.39667e+01,'min': 0, 'max':10, 'fixed':'d b0 2'}
conf['params']['collinspi']['s b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b1 2'}
conf['params']['collinspi']['s c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c0 2'}
conf['params']['collinspi']['s c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c1 2'}
conf['params']['collinspi']['s d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d0 2'}
conf['params']['collinspi']['s d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d1 2'}


conf['params']['collinspi']['ub N0 1'] ={'value':   -5.94721e+00,'min': -15, 'max': 10,  'fixed':'d N0 1'}
conf['params']['collinspi']['ub N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['collinspi']['ub a0 1'] ={'value':    4.25846e+00,'min':-1, 'max': 20,  'fixed':'d a0 1'}
conf['params']['collinspi']['ub a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['collinspi']['ub b0 1'] ={'value':    4.19000e+00,'min': 2.19, 'max':10,    'fixed':'d b0 1'}
conf['params']['collinspi']['ub b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 1'}
conf['params']['collinspi']['ub c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 1'}
conf['params']['collinspi']['ub c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 1'}
conf['params']['collinspi']['ub d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 1'}
conf['params']['collinspi']['ub d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 1'}
conf['params']['collinspi']['ub N0 2'] ={'value':   -3.06065e+02,'min': -10, 'max': 10, 'fixed':'d N0 2'}
conf['params']['collinspi']['ub N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'d N1 2'}
conf['params']['collinspi']['ub a0 2'] ={'value':    8.69912e-01,'min':-2, 'max': 20,   'fixed':'d a0 2'}
conf['params']['collinspi']['ub a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'d a1 2'}
conf['params']['collinspi']['ub b0 2'] ={'value':    1.39667e+01,'min': 0, 'max':10,    'fixed':'d b0 2'}
conf['params']['collinspi']['ub b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 2'}
conf['params']['collinspi']['ub c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 2'}
conf['params']['collinspi']['ub c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 2'}
conf['params']['collinspi']['ub d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 2'}
conf['params']['collinspi']['ub d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 2'}

conf['params']['collinspi']['db N0 1'] ={'value':    1.56490e+00,'min': -15, 'max': 15,  'fixed':'u N0 1'}
conf['params']['collinspi']['db N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'u N1 1'}
conf['params']['collinspi']['db a0 1'] ={'value':    1.23346e-01,'min':-0.228, 'max': 5,  'fixed':'u a0 1'}
conf['params']['collinspi']['db a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'u a1 1'}
conf['params']['collinspi']['db b0 1'] ={'value':    2.90807e+00,'min': 0.9, 'max':15,    'fixed':'u b0 1'}
conf['params']['collinspi']['db b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 1'}
conf['params']['collinspi']['db c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 1'}
conf['params']['collinspi']['db c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 1'}
conf['params']['collinspi']['db d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 1'}
conf['params']['collinspi']['db d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 1'}
conf['params']['collinspi']['db N0 2'] ={'value':    2.93621e+00,'min': -15, 'max': 15, 'fixed':'u N0 2'}
conf['params']['collinspi']['db N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'u N1 2'}
conf['params']['collinspi']['db a0 2'] ={'value':    7.95023e-01,'min':-0.228, 'max': 20,   'fixed':'u a0 2'}
conf['params']['collinspi']['db a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'u a1 2'}
conf['params']['collinspi']['db b0 2'] ={'value':    5.72097e+00,'min': 1.20, 'max':10,    'fixed':'u b0 2'}
conf['params']['collinspi']['db b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 2'}
conf['params']['collinspi']['db c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 2'}
conf['params']['collinspi']['db c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 2'}
conf['params']['collinspi']['db d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 2'}
conf['params']['collinspi']['db d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 2'}

conf['params']['collinspi']['sb N0 1'] ={'value':   -5.94721e+00,'min': -15, 'max': 10,  'fixed':'d N0 1'}
conf['params']['collinspi']['sb N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['collinspi']['sb a0 1'] ={'value':    4.25846e+00,'min':-1, 'max': 20,  'fixed':'d a0 1'}
conf['params']['collinspi']['sb a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['collinspi']['sb b0 1'] ={'value':    4.19000e+00,'min': 2.19, 'max':10,    'fixed':'d b0 1'}
conf['params']['collinspi']['sb b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 1'}
conf['params']['collinspi']['sb c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 1'}
conf['params']['collinspi']['sb c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 1'}
conf['params']['collinspi']['sb d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 1'}
conf['params']['collinspi']['sb d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 1'}
conf['params']['collinspi']['sb N0 2'] ={'value':   -3.06065e+02,'min': -10, 'max': 10, 'fixed':'d N0 2'}
conf['params']['collinspi']['sb N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'d N1 2'}
conf['params']['collinspi']['sb a0 2'] ={'value':    8.69912e-01,'min':-2, 'max': 20,   'fixed':'d a0 2'}
conf['params']['collinspi']['sb a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'d a1 2'}
conf['params']['collinspi']['sb b0 2'] ={'value':    1.39667e+01,'min': 0, 'max':10,    'fixed':'d b0 2'}
conf['params']['collinspi']['sb b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 2'}
conf['params']['collinspi']['sb c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 2'}
conf['params']['collinspi']['sb c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 2'}
conf['params']['collinspi']['sb d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 2'}
conf['params']['collinspi']['sb d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 2'}


#--steps
conf['steps']={}


##--unpol + cos2phi + sia Collins
conf['steps'][1]={}
conf['steps'][1]['dep']=[] #[1,6] #[2,4]
conf['steps'][1]['active distributions']=['pdf','ffpi','collinspi', 'boermulders','ffk']
conf['steps'][1]['datasets']={}

conf['steps'][1]['datasets']['sidis']=[]
conf['steps'][1]['datasets']['sidis'].append(1000) #'sidis/expdata/1000.xlsx'  # |  proton   | pi+   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1001) #'sidis/expdata/1001.xlsx'  # |  proton   | pi-   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1002) #'sidis/expdata/1002.xlsx'  # |  proton   | k+   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1003) #'sidis/expdata/1003.xlsx'  # |  proton   | k-   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1004) #'sidis/expdata/1004.xlsx'  # |  deuteron | pi+   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1005) #'sidis/expdata/1005.xlsx'  # |  deuteron | pi-   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1006) #'sidis/expdata/1004.xlsx'  # |  deuteron | k+   | M_Hermes | HERMES
conf['steps'][1]['datasets']['sidis'].append(1007) #'sidis/expdata/1005.xlsx'  # |  deuteron | k-   | M_Hermes | HERMES

conf['steps'][1]['datasets']['sidis'].append(5002) #'sidis/expdata/5002.xlsx'  # |  proton   | pi-   | AUUcos2 | HERMES
conf['steps'][1]['datasets']['sidis'].append(5003) #'sidis/expdata/5003.xlsx'  # |  proton   | pi+   | AUUcos2 | HERMES
conf['steps'][1]['datasets']['sidis'].append(5006) #'sidis/expdata/5006.xlsx'  # |  deuteron | pi-   | AUUcos2 | HERMES
conf['steps'][1]['datasets']['sidis'].append(5007) #'sidis/expdata/5007.xlsx'  # |  deuteron | pi+   | AUUcos2 | HERMES


conf['steps'][1]['datasets']['sia']=[]
conf['steps'][1]['datasets']['sia'].append(1000) # babar | pi,pi | AUL-0     | 9      | z1,z2,pT0  |
conf['steps'][1]['datasets']['sia'].append(1001) # babar | pi,pi | AUC-0     | 9      | z1,z2,pT0  |
conf['steps'][1]['datasets']['sia'].append(1002) # babar | pi,pi | AUC-0     | 36     | z1,z2      |
conf['steps'][1]['datasets']['sia'].append(1003) # babar | pi,pi | AUL-0     | 36     | z1,z2      |
conf['steps'][1]['datasets']['sia'].append(1004) # belle | pi,pi | AUT-0-CCP | 16     | z1,z2,qT   |
conf['steps'][1]['datasets']['sia'].append(1005) # belle | pi,pi | AUT-0     | 16     | z1,z2,qT   |
conf['steps'][1]['datasets']['sia'].append(2008) # babar | pi,pi | AUL-0     | 16     | z1,z2      |
conf['steps'][1]['datasets']['sia'].append(2009) # babar | pi,pi | AUC-0     | 16     | z1,z2      |
conf['steps'][1]['datasets']['sia'].append(3000) # bes3 | pi,pi | AUL-0     | 6     | z1,z2        |
conf['steps'][1]['datasets']['sia'].append(3001) # bes3 | pi,pi | AUC-0     | 6     | z1,z2        |
conf['steps'][1]['datasets']['sia'].append(3002) # bes3 | pi,pi | AUL-0     | 5     | z1,z2,pT     |
conf['steps'][1]['datasets']['sia'].append(3003) # bes3 | pi,pi | AUC-0     | 5     | z1,z2.pT     |
