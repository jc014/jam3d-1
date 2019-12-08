conf={}


############################################################################
# resouce allocation
conf['ncpus']=4

############################################################################
# maxlike setup
conf['screen mode']='plain'
#conf['screen mode']='curses'

#conf['shape'] = 1

# params
conf['params']={}
# Parameters in gaussian approximation, parton model:
# TMD PDF:
conf['params']['pdf']={}
conf['params']['pdf']['widths1_uv']  ={'value':    5.24140e-01,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths2_uv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['pdf']['widths1_dv']  ={'value':    5.24140e-01,'min': 0,'max':1,'fixed':'widths1_uv'}
conf['params']['pdf']['widths2_dv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}
conf['params']['pdf']['widths1_sea'] ={'value':    5.84650e-01,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths2_sea'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}

# TMD FF:
conf['params']['ffpi']={}
conf['params']['ffpi']['widths1_fav']   ={'value':    1.24050e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffpi']['widths1_ufav'] ={'value':    1.43730e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths2_ufav'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}
conf['params']['ffk']={}
conf['params']['ffk']['widths1_fav']    ={'value':    1.32333e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_fav']    ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']['widths1_ufav']  ={'value':    2.02660e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}
conf['params']['ffh']={}
conf['params']['ffh']['widths1_fav']    ={'value':    1.24050e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_fav']    ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']['widths1_ufav']  ={'value':    1.43730e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}


# Transversity
conf['params']['transversity']={}
conf['params']['transversity']['widths1_uv']  ={'value':    9.12760e-01, 'min':0, 'max':1, 'fixed': False}
conf['params']['transversity']['widths1_dv']  ={'value':    9.12760e-01, 'min':0, 'max':1, 'fixed': 'widths1_uv'}
conf['params']['transversity']['widths1_sea'] ={'value':    5.35004e-01, 'min':0, 'max':1, 'fixed': True}
conf['params']['transversity']['widths2_uv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': True}
conf['params']['transversity']['widths2_dv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}
conf['params']['transversity']['widths2_sea'] ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}

conf['params']['transversity']['u N0 1'] ={'value':    9.43653e+00, 'min': -10.0, 'max': 10.0, 'fixed': False}
conf['params']['transversity']['u N1 1'] ={'value':    0.00000e+00, 'min': -10.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u a0 1'] ={'value':    7.17307e-01, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['u a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['u b0 1'] ={'value':    2.54075e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['u b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['u c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['u a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['u b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['u b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['u c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['u d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['transversity']['d N0 1'] ={'value':   -2.27102e+01, 'min': -100.0, 'max': 100.0, 'fixed': False}
conf['params']['transversity']['d N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 'u N1 1'}
conf['params']['transversity']['d a0 1'] ={'value':    1.61348e+00, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['d a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 'u a1 1'}
conf['params']['transversity']['d b0 1'] ={'value':    4.23995e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['d b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 'u b1 1'}
conf['params']['transversity']['d c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['d a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['d b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['d b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['d c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['d d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}

conf['params']['transversity']['s N0 1'] ={'value':    0.00000e+00, 'min': -10.0, 'max': 10.0, 'fixed': True}
conf['params']['transversity']['s N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['transversity']['s c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['transversity']['s d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['transversity']['ub N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['transversity']['ub N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['transversity']['ub a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['transversity']['ub a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['transversity']['ub b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['transversity']['ub b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['transversity']['ub c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['transversity']['ub c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['transversity']['ub d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['transversity']['ub d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['transversity']['ub N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['transversity']['ub N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['transversity']['ub a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['transversity']['ub a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['transversity']['ub b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['transversity']['ub b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['transversity']['ub c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['transversity']['ub c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['transversity']['ub d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['transversity']['ub d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

conf['params']['transversity']['db N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['transversity']['db N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['transversity']['db a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['transversity']['db a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['transversity']['db b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['transversity']['db b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['transversity']['db c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['transversity']['db c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['transversity']['db d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['transversity']['db d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['transversity']['db N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['transversity']['db N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['transversity']['db a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['transversity']['db a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['transversity']['db b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['transversity']['db b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['transversity']['db c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['transversity']['db c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['transversity']['db d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['transversity']['db d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

conf['params']['transversity']['sb N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['transversity']['sb N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['transversity']['sb a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['transversity']['sb a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['transversity']['sb b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['transversity']['sb b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['transversity']['sb c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['transversity']['sb c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['transversity']['sb d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['transversity']['sb d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['transversity']['sb N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['transversity']['sb N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['transversity']['sb a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['transversity']['sb a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['transversity']['sb b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['transversity']['sb b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['transversity']['sb c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['transversity']['sb c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['transversity']['sb d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['transversity']['sb d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}



# Collins TMD FF:
conf['params']['collinspi']={}

conf['params']['collinspi']['widths1_fav']   ={'value':    4.33929e-02,'min':0, 'max':1, 'fixed':False}
conf['params']['collinspi']['widths1_ufav'] ={'value':    1.96171e-01,'min':0, 'max':1, 'fixed':False}
conf['params']['collinspi']['widths2_fav']   ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['collinspi']['widths2_ufav'] ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':'widths2_fav'}

conf['params']['collinspi']['u N0 1'] ={'value':    2.51599e+00,'min': -15, 'max': 10, 'fixed':False}
conf['params']['collinspi']['u N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':True}
conf['params']['collinspi']['u a0 1'] ={'value':   -8.00000e-01,'min':-0.8, 'max': 5, 'fixed':False}
conf['params']['collinspi']['u a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5, 'fixed':True}
conf['params']['collinspi']['u b0 1'] ={'value':    7.00924e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['u b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['collinspi']['u c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u N0 2'] ={'value':    1.33319e+01,'min': -15, 'max': 15, 'fixed':False}
conf['params']['collinspi']['u N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u a0 2'] ={'value':    3.89137e+00,'min':-2, 'max': 20, 'fixed':False}
conf['params']['collinspi']['u a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['collinspi']['u b0 2'] ={'value':    3.29884e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['u b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['collinspi']['u c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['collinspi']['d N0 1'] ={'value':   -4.08898e+00,'min': -10, 'max': 10, 'fixed':False}
conf['params']['collinspi']['d N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'u N1 1'}
conf['params']['collinspi']['d a0 1'] ={'value':   -8.00000e-01,'min':-0.8, 'max': 10, 'fixed':False}
conf['params']['collinspi']['d a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'u a1 1'}
conf['params']['collinspi']['d b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['d b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'u b1 1'}
conf['params']['collinspi']['d c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d N0 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['collinspi']['d a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['collinspi']['d b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['collinspi']['d b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['collinspi']['d c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['d d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['collinspi']['s N0 1'] ={'value':   -4.08898e+00,'min': 0, 'max': 1, 'fixed':'d N0 1'}
conf['params']['collinspi']['s N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 1'}
conf['params']['collinspi']['s a0 1'] ={'value':   -8.00000e-01,'min':-2, 'max': 2, 'fixed':'d a0 1'}
conf['params']['collinspi']['s a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 1'}
conf['params']['collinspi']['s b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10, 'fixed':'d b0 1'}
conf['params']['collinspi']['s b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b1 1'}
conf['params']['collinspi']['s c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c0 1'}
conf['params']['collinspi']['s c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c1 1'}
conf['params']['collinspi']['s d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d0 1'}
conf['params']['collinspi']['s d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d1 1'}
conf['params']['collinspi']['s N0 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N0 2'}
conf['params']['collinspi']['s N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 2'}
conf['params']['collinspi']['s a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a0 2'}
conf['params']['collinspi']['s a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 2'}
conf['params']['collinspi']['s b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b0 2'}
conf['params']['collinspi']['s b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b1 2'}
conf['params']['collinspi']['s c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c0 2'}
conf['params']['collinspi']['s c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c1 2'}
conf['params']['collinspi']['s d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d0 2'}
conf['params']['collinspi']['s d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d1 2'}


conf['params']['collinspi']['ub N0 1'] ={'value':   -4.08898e+00,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['collinspi']['ub N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['collinspi']['ub a0 1'] ={'value':   -8.00000e-01,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['collinspi']['ub a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['collinspi']['ub b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10,    'fixed':'d b0 1'}
conf['params']['collinspi']['ub b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 1'}
conf['params']['collinspi']['ub c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 1'}
conf['params']['collinspi']['ub c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 1'}
conf['params']['collinspi']['ub d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 1'}
conf['params']['collinspi']['ub d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 1'}
conf['params']['collinspi']['ub N0 2'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N0 2'}
conf['params']['collinspi']['ub N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'d N1 2'}
conf['params']['collinspi']['ub a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 20,   'fixed':'d a0 2'}
conf['params']['collinspi']['ub a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'d a1 2'}
conf['params']['collinspi']['ub b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b0 2'}
conf['params']['collinspi']['ub b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 2'}
conf['params']['collinspi']['ub c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 2'}
conf['params']['collinspi']['ub c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 2'}
conf['params']['collinspi']['ub d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 2'}
conf['params']['collinspi']['ub d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 2'}

conf['params']['collinspi']['db N0 1'] ={'value':    2.51599e+00,'min': -15, 'max': 2,  'fixed':'u N0 1'}
conf['params']['collinspi']['db N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'u N1 1'}
conf['params']['collinspi']['db a0 1'] ={'value':   -8.00000e-01,'min':-2.5, 'max': 5,  'fixed':'u a0 1'}
conf['params']['collinspi']['db a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'u a1 1'}
conf['params']['collinspi']['db b0 1'] ={'value':    7.00924e+00,'min': 0, 'max':10,    'fixed':'u b0 1'}
conf['params']['collinspi']['db b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 1'}
conf['params']['collinspi']['db c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 1'}
conf['params']['collinspi']['db c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 1'}
conf['params']['collinspi']['db d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 1'}
conf['params']['collinspi']['db d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 1'}
conf['params']['collinspi']['db N0 2'] ={'value':    1.33319e+01,'min': -10, 'max': 10, 'fixed':'u N0 2'}
conf['params']['collinspi']['db N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'u N1 2'}
conf['params']['collinspi']['db a0 2'] ={'value':    3.89137e+00,'min':-2, 'max': 20,   'fixed':'u a0 2'}
conf['params']['collinspi']['db a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'u a1 2'}
conf['params']['collinspi']['db b0 2'] ={'value':    3.29884e+00,'min': 0, 'max':10,    'fixed':'u b0 2'}
conf['params']['collinspi']['db b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 2'}
conf['params']['collinspi']['db c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 2'}
conf['params']['collinspi']['db c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 2'}
conf['params']['collinspi']['db d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 2'}
conf['params']['collinspi']['db d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 2'}

conf['params']['collinspi']['sb N0 1'] ={'value':   -4.08898e+00,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['collinspi']['sb N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['collinspi']['sb a0 1'] ={'value':   -8.00000e-01,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['collinspi']['sb a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['collinspi']['sb b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10,    'fixed':'d b0 1'}
conf['params']['collinspi']['sb b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 1'}
conf['params']['collinspi']['sb c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 1'}
conf['params']['collinspi']['sb c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 1'}
conf['params']['collinspi']['sb d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 1'}
conf['params']['collinspi']['sb d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 1'}
conf['params']['collinspi']['sb N0 2'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N0 2'}
conf['params']['collinspi']['sb N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'d N1 2'}
conf['params']['collinspi']['sb a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 20,   'fixed':'d a0 2'}
conf['params']['collinspi']['sb a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'d a1 2'}
conf['params']['collinspi']['sb b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b0 2'}
conf['params']['collinspi']['sb b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 2'}
conf['params']['collinspi']['sb c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 2'}
conf['params']['collinspi']['sb c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 2'}
conf['params']['collinspi']['sb d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 2'}
conf['params']['collinspi']['sb d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 2'}

#Htilde FF
conf['params']['Htildepi']={}
conf['params']['Htildepi']['widths1_fav']   ={'value':    1.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['Htildepi']['widths1_ufav'] ={'value':    1.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['Htildepi']['widths2_fav']   ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['Htildepi']['widths2_ufav'] ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}

conf['params']['Htildepi']['u N0 1'] ={'value':   -2.74176e-01,'min': -15, 'max': 10, 'fixed':False}
conf['params']['Htildepi']['u N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':True}
conf['params']['Htildepi']['u a0 1'] ={'value':    4.29342e+00,'min':-2, 'max': 5, 'fixed':False}
conf['params']['Htildepi']['u a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5, 'fixed':True}
conf['params']['Htildepi']['u b0 1'] ={'value':    1.65647e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['Htildepi']['u b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['Htildepi']['u c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u N0 2'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':True}
conf['params']['Htildepi']['u N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 20, 'fixed':True}
conf['params']['Htildepi']['u a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['Htildepi']['u b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['Htildepi']['u b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['Htildepi']['u c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['u d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['Htildepi']['d N0 1'] ={'value':    1.52788e+00,'min': -10, 'max': 15, 'fixed':False}
conf['params']['Htildepi']['d N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'u N1 1'}
conf['params']['Htildepi']['d a0 1'] ={'value':    3.73506e+00,'min':-2, 'max': 5, 'fixed':False}
conf['params']['Htildepi']['d a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'u a1 1'}
conf['params']['Htildepi']['d b0 1'] ={'value':    1.83252e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['Htildepi']['d b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'u b1 1'}
conf['params']['Htildepi']['d c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d N0 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['Htildepi']['d a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['Htildepi']['d b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['Htildepi']['d b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
conf['params']['Htildepi']['d c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['Htildepi']['s N0 1'] ={'value':    1.52788e+00,'min': 0, 'max': 1, 'fixed':'d N0 1'}
conf['params']['Htildepi']['s N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 1'}
conf['params']['Htildepi']['s a0 1'] ={'value':    3.73506e+00,'min':-2, 'max': 2, 'fixed':'d a0 1'}
conf['params']['Htildepi']['s a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 1'}
conf['params']['Htildepi']['s b0 1'] ={'value':    1.83252e+00,'min': 0, 'max':10, 'fixed':'d b0 1'}
conf['params']['Htildepi']['s b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b1 1'}
conf['params']['Htildepi']['s c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c0 1'}
conf['params']['Htildepi']['s c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c1 1'}
conf['params']['Htildepi']['s d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d0 1'}
conf['params']['Htildepi']['s d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d1 1'}
conf['params']['Htildepi']['s N0 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N0 2'}
conf['params']['Htildepi']['s N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 2'}
conf['params']['Htildepi']['s a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a0 2'}
conf['params']['Htildepi']['s a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 2'}
conf['params']['Htildepi']['s b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b0 2'}
conf['params']['Htildepi']['s b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':'d b1 2'}
conf['params']['Htildepi']['s c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c0 2'}
conf['params']['Htildepi']['s c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d c1 2'}
conf['params']['Htildepi']['s d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d0 2'}
conf['params']['Htildepi']['s d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':'d d1 2'}


conf['params']['Htildepi']['ub N0 1'] ={'value':    1.52788e+00,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['Htildepi']['ub N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['Htildepi']['ub a0 1'] ={'value':    3.73506e+00,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['Htildepi']['ub a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['Htildepi']['ub b0 1'] ={'value':    1.83252e+00,'min': 0, 'max':10,    'fixed':'d b0 1'}
conf['params']['Htildepi']['ub b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 1'}
conf['params']['Htildepi']['ub c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 1'}
conf['params']['Htildepi']['ub c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 1'}
conf['params']['Htildepi']['ub d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 1'}
conf['params']['Htildepi']['ub d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 1'}
conf['params']['Htildepi']['ub N0 2'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N0 2'}
conf['params']['Htildepi']['ub N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'d N1 2'}
conf['params']['Htildepi']['ub a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 20,   'fixed':'d a0 2'}
conf['params']['Htildepi']['ub a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'d a1 2'}
conf['params']['Htildepi']['ub b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b0 2'}
conf['params']['Htildepi']['ub b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 2'}
conf['params']['Htildepi']['ub c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 2'}
conf['params']['Htildepi']['ub c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 2'}
conf['params']['Htildepi']['ub d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 2'}
conf['params']['Htildepi']['ub d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 2'}

conf['params']['Htildepi']['db N0 1'] ={'value':   -2.74176e-01,'min': -15, 'max': 2,  'fixed':'u N0 1'}
conf['params']['Htildepi']['db N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'u N1 1'}
conf['params']['Htildepi']['db a0 1'] ={'value':    4.29342e+00,'min':-2.5, 'max': 5,  'fixed':'u a0 1'}
conf['params']['Htildepi']['db a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'u a1 1'}
conf['params']['Htildepi']['db b0 1'] ={'value':    1.65647e+00,'min': 0, 'max':10,    'fixed':'u b0 1'}
conf['params']['Htildepi']['db b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 1'}
conf['params']['Htildepi']['db c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 1'}
conf['params']['Htildepi']['db c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 1'}
conf['params']['Htildepi']['db d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 1'}
conf['params']['Htildepi']['db d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 1'}
conf['params']['Htildepi']['db N0 2'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'u N0 2'}
conf['params']['Htildepi']['db N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'u N1 2'}
conf['params']['Htildepi']['db a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 20,   'fixed':'u a0 2'}
conf['params']['Htildepi']['db a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'u a1 2'}
conf['params']['Htildepi']['db b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b0 2'}
conf['params']['Htildepi']['db b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 2'}
conf['params']['Htildepi']['db c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 2'}
conf['params']['Htildepi']['db c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 2'}
conf['params']['Htildepi']['db d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 2'}
conf['params']['Htildepi']['db d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 2'}

conf['params']['Htildepi']['sb N0 1'] ={'value':    1.52788e+00,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['Htildepi']['sb N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['Htildepi']['sb a0 1'] ={'value':    3.73506e+00,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['Htildepi']['sb a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['Htildepi']['sb b0 1'] ={'value':    1.83252e+00,'min': 0, 'max':10,    'fixed':'d b0 1'}
conf['params']['Htildepi']['sb b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 1'}
conf['params']['Htildepi']['sb c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 1'}
conf['params']['Htildepi']['sb c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 1'}
conf['params']['Htildepi']['sb d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 1'}
conf['params']['Htildepi']['sb d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 1'}
conf['params']['Htildepi']['sb N0 2'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N0 2'}
conf['params']['Htildepi']['sb N1 2'] ={'value':    0.00000e+00,'min': 0, 'max': 1,    'fixed':'d N1 2'}
conf['params']['Htildepi']['sb a0 2'] ={'value':    0.00000e+00,'min':-2, 'max': 20,   'fixed':'d a0 2'}
conf['params']['Htildepi']['sb a1 2'] ={'value':    0.00000e+00,'min':-2, 'max': 2,    'fixed':'d a1 2'}
conf['params']['Htildepi']['sb b0 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b0 2'}
conf['params']['Htildepi']['sb b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'d b1 2'}
conf['params']['Htildepi']['sb c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c0 2'}
conf['params']['Htildepi']['sb c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d c1 2'}
conf['params']['Htildepi']['sb d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d0 2'}
conf['params']['Htildepi']['sb d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'d d1 2'}

#Sivers function
conf['params']['sivers']={}
conf['params']['sivers']['widths1_uv']  ={'value':    2.01072e-01, 'min':0, 'max':1, 'fixed': False}
conf['params']['sivers']['widths1_dv']  ={'value':    2.01072e-01, 'min':0, 'max':1, 'fixed': 'widths1_uv'}
conf['params']['sivers']['widths1_sea'] ={'value':    3.88407e-01, 'min':0, 'max':1, 'fixed': True}
conf['params']['sivers']['widths2_uv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': True}
conf['params']['sivers']['widths2_dv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}
conf['params']['sivers']['widths2_sea'] ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}

conf['params']['sivers']['u N0 1'] ={'value':   -1.82741e-01, 'min': -10.0, 'max': 10.0, 'fixed': False}
conf['params']['sivers']['u N1 1'] ={'value':    0.00000e+00, 'min': -10.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u a0 1'] ={'value':   -9.26734e-02, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['sivers']['u a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['u b0 1'] ={'value':    4.42200e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
conf['params']['sivers']['u b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['u c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['u a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['u b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['u b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['u c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['u d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['sivers']['d N0 1'] ={'value':    5.85923e-02, 'min': -20.0, 'max': 10.0, 'fixed': False}
conf['params']['sivers']['d N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 'u N1 1'}
conf['params']['sivers']['d a0 1'] ={'value':   -5.87940e-01, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['sivers']['d a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 'u a1 1'}
conf['params']['sivers']['d b0 1'] ={'value':    3.41349e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
conf['params']['sivers']['d b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 'u b1 1'}
conf['params']['sivers']['d c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['d a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['d b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['d b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['d c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['d d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['sivers']['s N0 1'] ={'value':    0.00000e+00, 'min': -10.0, 'max': 10.0, 'fixed': True}
conf['params']['sivers']['s N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': True}
conf['params']['sivers']['s c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}
conf['params']['sivers']['s d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': True}


conf['params']['sivers']['ub N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['sivers']['ub N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['sivers']['ub a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['sivers']['ub a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['sivers']['ub b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['sivers']['ub b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['sivers']['ub c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['sivers']['ub c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['sivers']['ub d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['sivers']['ub d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['sivers']['ub N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['sivers']['ub N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['sivers']['ub a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['sivers']['ub a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['sivers']['ub b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['sivers']['ub b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['sivers']['ub c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['sivers']['ub c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['sivers']['ub d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['sivers']['ub d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

conf['params']['sivers']['db N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['sivers']['db N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['sivers']['db a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['sivers']['db a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['sivers']['db b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['sivers']['db b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['sivers']['db c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['sivers']['db c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['sivers']['db d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['sivers']['db d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['sivers']['db N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['sivers']['db N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['sivers']['db a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['sivers']['db a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['sivers']['db b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['sivers']['db b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['sivers']['db c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['sivers']['db c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['sivers']['db d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['sivers']['db d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

conf['params']['sivers']['sb N0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 10 , 'fixed': 's N0 1'}
conf['params']['sivers']['sb N1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 1'}
conf['params']['sivers']['sb a0 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 1'}
conf['params']['sivers']['sb a1 1'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 1'}
conf['params']['sivers']['sb b0 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 1'}
conf['params']['sivers']['sb b1 1'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 1'}
conf['params']['sivers']['sb c0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 1'}
conf['params']['sivers']['sb c1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 1'}
conf['params']['sivers']['sb d0 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 1'}
conf['params']['sivers']['sb d1 1'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 1'}
conf['params']['sivers']['sb N0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N0 2'}
conf['params']['sivers']['sb N1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max': 1.0, 'fixed': 's N1 2'}
conf['params']['sivers']['sb a0 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a0 2'}
conf['params']['sivers']['sb a1 2'] ={'value':    0.00000e+00, 'min':-2.0, 'max':10.0, 'fixed': 's a1 2'}
conf['params']['sivers']['sb b0 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b0 2'}
conf['params']['sivers']['sb b1 2'] ={'value':    0.00000e+00, 'min': 0.0, 'max':10.0, 'fixed': 's b1 2'}
conf['params']['sivers']['sb c0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c0 2'}
conf['params']['sivers']['sb c1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's c1 2'}
conf['params']['sivers']['sb d0 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d0 2'}
conf['params']['sivers']['sb d1 2'] ={'value':    0.00000e+00, 'min':-1.0, 'max': 1.0, 'fixed': 's d1 2'}

############################################################################
# set data sets

conf['datasets']={}

# AN

conf['datasets']['ANep']={}


conf['datasets']['ANep']['norm']={}

conf['datasets']['ANep']['filters']=[]


conf['datasets']['ANep']['xlsx']={}

conf['datasets']['ANep']={}
conf['datasets']['ANep']['norm']={}
conf['datasets']['ANep']['xlsx']={}
conf['datasets']['ANep']['filters']={}

conf['datasets']['ANep']['filters']=['pT>1.0']

conf['datasets']['ANep']['xlsx'][1000]='AN_ep/expdata/1000.xlsx' # pi- AN-xFbins  pT
conf['datasets']['ANep']['xlsx'][1001]='AN_ep/expdata/1001.xlsx' # pi- AN-pTbins  xF
conf['datasets']['ANep']['xlsx'][2000]='AN_ep/expdata/2000.xlsx' # pi+ AN-xFbins  pT
conf['datasets']['ANep']['xlsx'][2001]='AN_ep/expdata/2001.xlsx' # pi+ AN-pTbins  xF

for k in conf['datasets']['ANep']['xlsx']: conf['datasets']['ANep']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}
