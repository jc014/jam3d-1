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
conf['params']['pdf']['widths2_sea'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}

# TMD FF:
conf['params']['ffpi']={}
conf['params']['ffpi']['widths1_fav']   ={'value':    1.24050e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffpi']['widths1_ufav'] ={'value':    1.43730e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths2_ufav'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']={}
conf['params']['ffk']['widths1_fav']    ={'value':    1.32333e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_fav']    ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']['widths1_ufav']  ={'value':    2.02660e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']={}
conf['params']['ffh']['widths1_fav']    ={'value':    1.24050e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_fav']    ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']['widths1_ufav']  ={'value':    1.43730e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}


# Transversity
conf['params']['transversity']={}
conf['params']['transversity']['widths1_uv']  ={'value':    6.36914e-01, 'min':0, 'max':1, 'fixed': False}
conf['params']['transversity']['widths1_dv']  ={'value':    6.36914e-01, 'min':0, 'max':1, 'fixed': 'widths1_uv'}
conf['params']['transversity']['widths1_sea'] ={'value':    5.35004e-01, 'min':0, 'max':1, 'fixed': False}
conf['params']['transversity']['widths2_uv']  ={'value':    1.50864e-02, 'min':0, 'max':1, 'fixed': False}
conf['params']['transversity']['widths2_dv']  ={'value':    1.50864e-02, 'min':0, 'max':1, 'fixed': 'widths2_uv'}
conf['params']['transversity']['widths2_sea'] ={'value':    1.50864e-02, 'min':0, 'max':1, 'fixed': 'widths2_uv'}

conf['params']['transversity']['u N0 1'] ={'value':    9.24040e+00, 'min': -10.0, 'max': 10.0, 'fixed': False}
conf['params']['transversity']['u N1 1'] ={'value':   -4.75538e-01, 'min': -10.0, 'max': 1.0, 'fixed': False}
conf['params']['transversity']['u a0 1'] ={'value':    4.91144e-01, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['u a1 1'] ={'value':    7.25593e-02, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['u b0 1'] ={'value':    2.72459e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['u b1 1'] ={'value':    1.31194e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
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


conf['params']['transversity']['d N0 1'] ={'value':   -2.28339e+01, 'min': -100.0, 'max': 100.0, 'fixed': False}
conf['params']['transversity']['d N1 1'] ={'value':   -4.75538e-01, 'min': 0.0, 'max': 1.0, 'fixed': 'u N1 1'}
conf['params']['transversity']['d a0 1'] ={'value':    1.68378e+00, 'min':-2.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['d a1 1'] ={'value':    7.25593e-02, 'min':-2.0, 'max':10.0, 'fixed': 'u a1 1'}
conf['params']['transversity']['d b0 1'] ={'value':    4.24694e+00, 'min': 0.0, 'max':10.0, 'fixed': False}
conf['params']['transversity']['d b1 1'] ={'value':    1.31194e+00, 'min': 0.0, 'max':10.0, 'fixed': 'u b1 1'}
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
conf['params']['collinspi']['widths2_fav']   ={'value':    2.88646e-08,'min':0, 'max':1, 'fixed':False}
conf['params']['collinspi']['widths2_ufav'] ={'value':    2.88646e-08,'min':0, 'max':1, 'fixed':'widths2_fav'}

conf['params']['collinspi']['u N0 1'] ={'value':    2.51599e+00,'min': -15, 'max': 10, 'fixed':False}
conf['params']['collinspi']['u N1 1'] ={'value':    1.28385e-07,'min': -10, 'max': 10, 'fixed':False}
conf['params']['collinspi']['u a0 1'] ={'value':   -8.00000e-01,'min':-0.8, 'max': 5, 'fixed':False}
conf['params']['collinspi']['u a1 1'] ={'value':   -4.66241e-08,'min':-2, 'max': 5, 'fixed':False}
conf['params']['collinspi']['u b0 1'] ={'value':    7.00924e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['u b1 1'] ={'value':    3.64940e-08,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['u c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u N0 2'] ={'value':    1.33319e+01,'min': -15, 'max': 15, 'fixed':False}
conf['params']['collinspi']['u N1 2'] ={'value':    7.50293e-10,'min': 0, 'max': 1, 'fixed':False}
conf['params']['collinspi']['u a0 2'] ={'value':    3.89137e+00,'min':-2, 'max': 20, 'fixed':False}
conf['params']['collinspi']['u a1 2'] ={'value':   -1.20988e-09,'min':-2, 'max': 2, 'fixed':False}
conf['params']['collinspi']['u b0 2'] ={'value':    3.29884e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['u b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['u c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}
conf['params']['collinspi']['u d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1, 'fixed':True}

conf['params']['collinspi']['d N0 1'] ={'value':   -4.08898e+00,'min': -10, 'max': 10, 'fixed':False}
conf['params']['collinspi']['d N1 1'] ={'value':    1.28385e-07,'min': 0, 'max': 1, 'fixed':'u N1 1'}
conf['params']['collinspi']['d a0 1'] ={'value':   -8.00000e-01,'min':-0.8, 'max': 10, 'fixed':False}
conf['params']['collinspi']['d a1 1'] ={'value':   -4.66241e-08,'min':-2, 'max': 2, 'fixed':'u a1 1'}
conf['params']['collinspi']['d b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10, 'fixed':False}
conf['params']['collinspi']['d b1 1'] ={'value':    3.64940e-08,'min': 0, 'max':10, 'fixed':'u b1 1'}
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
conf['params']['collinspi']['s N1 1'] ={'value':    1.28385e-07,'min': 0, 'max': 1, 'fixed':'d N1 1'}
conf['params']['collinspi']['s a0 1'] ={'value':   -8.00000e-01,'min':-2, 'max': 2, 'fixed':'d a0 1'}
conf['params']['collinspi']['s a1 1'] ={'value':   -4.66241e-08,'min':-2, 'max': 2, 'fixed':'d a1 1'}
conf['params']['collinspi']['s b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10, 'fixed':'d b0 1'}
conf['params']['collinspi']['s b1 1'] ={'value':    3.64940e-08,'min': 0, 'max':10, 'fixed':'d b1 1'}
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
conf['params']['collinspi']['ub N1 1'] ={'value':    1.28385e-07,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['collinspi']['ub a0 1'] ={'value':   -8.00000e-01,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['collinspi']['ub a1 1'] ={'value':   -4.66241e-08,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['collinspi']['ub b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10,    'fixed':'d b0 1'}
conf['params']['collinspi']['ub b1 1'] ={'value':    3.64940e-08,'min': 0, 'max':10,    'fixed':'d b1 1'}
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
conf['params']['collinspi']['db N1 1'] ={'value':    1.28385e-07,'min': -10, 'max': 10, 'fixed':'u N1 1'}
conf['params']['collinspi']['db a0 1'] ={'value':   -8.00000e-01,'min':-2.5, 'max': 5,  'fixed':'u a0 1'}
conf['params']['collinspi']['db a1 1'] ={'value':   -4.66241e-08,'min':-2, 'max': 5,    'fixed':'u a1 1'}
conf['params']['collinspi']['db b0 1'] ={'value':    7.00924e+00,'min': 0, 'max':10,    'fixed':'u b0 1'}
conf['params']['collinspi']['db b1 1'] ={'value':    3.64940e-08,'min': 0, 'max':10,    'fixed':'u b1 1'}
conf['params']['collinspi']['db c0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 1'}
conf['params']['collinspi']['db c1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 1'}
conf['params']['collinspi']['db d0 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 1'}
conf['params']['collinspi']['db d1 1'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 1'}
conf['params']['collinspi']['db N0 2'] ={'value':    1.33319e+01,'min': -10, 'max': 10, 'fixed':'u N0 2'}
conf['params']['collinspi']['db N1 2'] ={'value':    7.50293e-10,'min': 0, 'max': 1,    'fixed':'u N1 2'}
conf['params']['collinspi']['db a0 2'] ={'value':    3.89137e+00,'min':-2, 'max': 20,   'fixed':'u a0 2'}
conf['params']['collinspi']['db a1 2'] ={'value':   -1.20988e-09,'min':-2, 'max': 2,    'fixed':'u a1 2'}
conf['params']['collinspi']['db b0 2'] ={'value':    3.29884e+00,'min': 0, 'max':10,    'fixed':'u b0 2'}
conf['params']['collinspi']['db b1 2'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b1 2'}
conf['params']['collinspi']['db c0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c0 2'}
conf['params']['collinspi']['db c1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u c1 2'}
conf['params']['collinspi']['db d0 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d0 2'}
conf['params']['collinspi']['db d1 2'] ={'value':    0.00000e+00,'min':-1, 'max': 1,    'fixed':'u d1 2'}

conf['params']['collinspi']['sb N0 1'] ={'value':   -4.08898e+00,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['collinspi']['sb N1 1'] ={'value':    1.28385e-07,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['collinspi']['sb a0 1'] ={'value':   -8.00000e-01,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['collinspi']['sb a1 1'] ={'value':   -4.66241e-08,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['collinspi']['sb b0 1'] ={'value':    6.16509e+00,'min': 0, 'max':10,    'fixed':'d b0 1'}
conf['params']['collinspi']['sb b1 1'] ={'value':    3.64940e-08,'min': 0, 'max':10,    'fixed':'d b1 1'}
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

#mimic Collins, widths not required, but only 1 shape
conf['params']['Htildepi']={}
conf['params']['Htildepi']['widths1_fav']   ={'value':    1.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['Htildepi']['widths1_ufav'] ={'value':    1.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['Htildepi']['widths2_fav']   ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['Htildepi']['widths2_ufav'] ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}

conf['params']['Htildepi']['u N0 1'] ={'value':    4.05868e-02,'min': -15, 'max': 10, 'fixed':True}
conf['params']['Htildepi']['u N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':True}
conf['params']['Htildepi']['u a0 1'] ={'value':    4.13081e-01,'min':-2, 'max': 5, 'fixed':True}
conf['params']['Htildepi']['u a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5, 'fixed':True}
conf['params']['Htildepi']['u b0 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
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

conf['params']['Htildepi']['d N0 1'] ={'value':    7.16070e-02,'min': -10, 'max': 10, 'fixed':True}
conf['params']['Htildepi']['d N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':True}
conf['params']['Htildepi']['d a0 1'] ={'value':   -1.13844e-01,'min':-2, 'max': 2, 'fixed':True}
conf['params']['Htildepi']['d a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':True}
conf['params']['Htildepi']['d b0 1'] ={'value':    4.91011e-03,'min': 0, 'max':10, 'fixed':True}
conf['params']['Htildepi']['d b1 1'] ={'value':    0.00000e+00,'min': 0, 'max':10, 'fixed':True}
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

conf['params']['Htildepi']['s N0 1'] ={'value':    7.16070e-02,'min': 0, 'max': 1, 'fixed':'d N0 1'}
conf['params']['Htildepi']['s N1 1'] ={'value':    0.00000e+00,'min': 0, 'max': 1, 'fixed':'d N1 1'}
conf['params']['Htildepi']['s a0 1'] ={'value':   -1.13844e-01,'min':-2, 'max': 2, 'fixed':'d a0 1'}
conf['params']['Htildepi']['s a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 2, 'fixed':'d a1 1'}
conf['params']['Htildepi']['s b0 1'] ={'value':    4.91011e-03,'min': 0, 'max':10, 'fixed':'d b0 1'}
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


conf['params']['Htildepi']['ub N0 1'] ={'value':    7.16070e-02,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['Htildepi']['ub N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['Htildepi']['ub a0 1'] ={'value':   -1.13844e-01,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['Htildepi']['ub a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['Htildepi']['ub b0 1'] ={'value':    4.91011e-03,'min': 0, 'max':10,    'fixed':'d b0 1'}
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

conf['params']['Htildepi']['db N0 1'] ={'value':    4.05868e-02,'min': -15, 'max': 2,  'fixed':'u N0 1'}
conf['params']['Htildepi']['db N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'u N1 1'}
conf['params']['Htildepi']['db a0 1'] ={'value':    4.13081e-01,'min':-2.5, 'max': 5,  'fixed':'u a0 1'}
conf['params']['Htildepi']['db a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'u a1 1'}
conf['params']['Htildepi']['db b0 1'] ={'value':    0.00000e+00,'min': 0, 'max':10,    'fixed':'u b0 1'}
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

conf['params']['Htildepi']['sb N0 1'] ={'value':    7.16070e-02,'min': -15, 'max': 2,  'fixed':'d N0 1'}
conf['params']['Htildepi']['sb N1 1'] ={'value':    0.00000e+00,'min': -10, 'max': 10, 'fixed':'d N1 1'}
conf['params']['Htildepi']['sb a0 1'] ={'value':   -1.13844e-01,'min':-2.5, 'max': 5,  'fixed':'d a0 1'}
conf['params']['Htildepi']['sb a1 1'] ={'value':    0.00000e+00,'min':-2, 'max': 5,    'fixed':'d a1 1'}
conf['params']['Htildepi']['sb b0 1'] ={'value':    4.91011e-03,'min': 0, 'max':10,    'fixed':'d b0 1'}
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


############################################################################
# set data sets

conf['datasets']={}

# SIA

conf['datasets']['sia']={}


conf['datasets']['sia']['norm']={}

conf['datasets']['sia']['filters']=[] # npts    = 122 chi2    = 69.799935
#conf['datasets']['sia']['filters'].append("z1<0.7")
#conf['datasets']['sia']['filters'].append("z2<0.7")
#conf['datasets']['sia']['filters'].append("Q2>1.69")
#conf['datasets']['sia']['filters'].append("pT>0.2 and pT<0.9")


conf['datasets']['sia']['xlsx']={}

conf['datasets']['sia']['xlsx'][1000]='sia/expdata/1000.xlsx' # babar | pi,pi | AUL-0     | 9      | z1,z2,pT0  |
conf['datasets']['sia']['xlsx'][1001]='sia/expdata/1001.xlsx' # babar | pi,pi | AUC-0     | 9      | z1,z2,pT0  |
conf['datasets']['sia']['xlsx'][1002]='sia/expdata/1002.xlsx' # babar | pi,pi | AUC-0     | 36     | z1,z2      |
conf['datasets']['sia']['xlsx'][1003]='sia/expdata/1003.xlsx' # babar | pi,pi | AUL-0     | 36     | z1,z2      |
conf['datasets']['sia']['xlsx'][1004]='sia/expdata/1004.xlsx' # belle | pi,pi | AUT-0-CCP | 16     | z1,z2,qT   |
conf['datasets']['sia']['xlsx'][1005]='sia/expdata/1005.xlsx' # belle | pi,pi | AUT-0     | 16     | z1,z2,qT   |
conf['datasets']['sia']['xlsx'][2008]='sia/expdata/2008.xlsx' # | babar      | pi,pi | AUL-0            | 16     | z1,z2      |
conf['datasets']['sia']['xlsx'][2009]='sia/expdata/2009.xlsx' # | babar      | pi,pi | AUC-0            | 16     | z1,z2      |

for k in conf['datasets']['sia']['xlsx']: conf['datasets']['sia']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}

# Collins Asy
conf['datasets']['sidis']={}
conf['datasets']['sidis']['norm']={}
conf['datasets']['sidis']['xlsx']={}
conf['datasets']['sidis']['filters']={}

conf['datasets']['sidis']['filters']=["z>0.2 and z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9"]

#multiplicity?

#conf["datasets"]["sidis"]["xlsx"][4007]="sidis/expdata/4007.xlsx"  #  compass  deuteron   k+   pT
#conf["datasets"]["sidis"]["xlsx"][4006]="sidis/expdata/4006.xlsx"  #  compass  deuteron   k+    x
#conf["datasets"]["sidis"]["xlsx"][4008]="sidis/expdata/4008.xlsx"  #  compass  deuteron   k+    z
#conf["datasets"]["sidis"]["xlsx"][4010]="sidis/expdata/4010.xlsx"  #  compass  deuteron   k-   pT
#conf["datasets"]["sidis"]["xlsx"][4009]="sidis/expdata/4009.xlsx"  #  compass  deuteron   k-    x
#conf["datasets"]["sidis"]["xlsx"][4011]="sidis/expdata/4011.xlsx"  #  compass  deuteron   k-    z
conf["datasets"]["sidis"]["xlsx"][4001]="sidis/expdata/4001.xlsx"  #  compass  deuteron  pi+  pT
conf["datasets"]["sidis"]["xlsx"][4000]="sidis/expdata/4000.xlsx"  #  compass  deuteron  pi+   x
conf["datasets"]["sidis"]["xlsx"][4002]="sidis/expdata/4002.xlsx"  #  compass  deuteron  pi+   z
conf["datasets"]["sidis"]["xlsx"][4004]="sidis/expdata/4004.xlsx"  #  compass  deuteron  pi-   pT
conf["datasets"]["sidis"]["xlsx"][4003]="sidis/expdata/4003.xlsx"  #  compass  deuteron  pi-    x
conf["datasets"]["sidis"]["xlsx"][4005]="sidis/expdata/4005.xlsx"  #  compass  deuteron  pi-    z
#conf["datasets"]["sidis"]["xlsx"][6003]="sidis/expdata/6003.xlsx"  #  compass    proton   k+   pt
#conf["datasets"]["sidis"]["xlsx"][6004]="sidis/expdata/6004.xlsx"  #  compass    proton   k+    x
#conf["datasets"]["sidis"]["xlsx"][6005]="sidis/expdata/6005.xlsx"  #  compass    proton   k+    z
#conf["datasets"]["sidis"]["xlsx"][6000]="sidis/expdata/6000.xlsx"  #  compass    proton   k-   pt
#conf["datasets"]["sidis"]["xlsx"][6001]="sidis/expdata/6001.xlsx"  #  compass    proton   k-    x
#conf["datasets"]["sidis"]["xlsx"][6002]="sidis/expdata/6002.xlsx"  #  compass    proton   k-    z
conf["datasets"]["sidis"]["xlsx"][3027]="sidis/expdata/3027.xlsx"  #  compass    proton    pt
conf["datasets"]["sidis"]["xlsx"][3025]="sidis/expdata/3025.xlsx"  #  compass    proton     x
conf["datasets"]["sidis"]["xlsx"][3010]="sidis/expdata/3010.xlsx"  #  compass    proton     z
conf["datasets"]["sidis"]["xlsx"][3012]="sidis/expdata/3012.xlsx"  #  compass    proton  pi-   pt
conf["datasets"]["sidis"]["xlsx"][3005]="sidis/expdata/3005.xlsx"  #  compass    proton  pi-    x
conf["datasets"]["sidis"]["xlsx"][3013]="sidis/expdata/3013.xlsx"  #  compass    proton  pi-    z
#conf["datasets"]["sidis"]["xlsx"][3024]="sidis/expdata/3024.xlsx"  #   HERMES    proton   k+   pt
#conf["datasets"]["sidis"]["xlsx"][3007]="sidis/expdata/3007.xlsx"  #   HERMES    proton   k+    x
#conf["datasets"]["sidis"]["xlsx"][3008]="sidis/expdata/3008.xlsx"  #   HERMES    proton   k+    z
#conf["datasets"]["sidis"]["xlsx"][3021]="sidis/expdata/3021.xlsx"  #   HERMES    proton   k-   pt
#conf["datasets"]["sidis"]["xlsx"][3017]="sidis/expdata/3017.xlsx"  #   HERMES    proton   k-    x
#conf["datasets"]["sidis"]["xlsx"][3023]="sidis/expdata/3023.xlsx"  #   HERMES    proton   k-    z
conf["datasets"]["sidis"]["xlsx"][3026]="sidis/expdata/3026.xlsx"  #   HERMES    proton pi+   pt
conf["datasets"]["sidis"]["xlsx"][3000]="sidis/expdata/3000.xlsx"  #   HERMES    proton pi+    x
conf["datasets"]["sidis"]["xlsx"][3003]="sidis/expdata/3003.xlsx"  #   HERMES    proton pi+    z
conf["datasets"]["sidis"]["xlsx"][3016]="sidis/expdata/3016.xlsx"  #   HERMES    proton  pi-   pt
conf["datasets"]["sidis"]["xlsx"][3004]="sidis/expdata/3004.xlsx"  #   HERMES    proton  pi-    x
conf["datasets"]["sidis"]["xlsx"][3018]="sidis/expdata/3018.xlsx"  #   HERMES    proton  pi-    z

#conf["datasets"]["sidis"]["xlsx"][9011]="sidis/expdata/9011.xlsx"  #   COMPASS    proton  h+    z
#conf["datasets"]["sidis"]["xlsx"][9022]="sidis/expdata/9022.xlsx"  #   COMPASS    proton  h-    z
#conf["datasets"]["sidis"]["xlsx"][9033]="sidis/expdata/9033.xlsx"  #   COMPASS    proton  h+    x
#conf["datasets"]["sidis"]["xlsx"][9044]="sidis/expdata/9044.xlsx"  #   COMPASS    proton  h-    x

for k in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}


# Lattice

#conf['datasets']['moments']={}
#conf['datasets']['moments']['filters']=[]
#conf['datasets']['moments']['filters'][0]['idx']=[1000,1001,1002,1002,1003,1004,1005,1006]

#conf['datasets']['moments']['xlsx']={}
#conf['datasets']['moments']['xlsx'][1000]='lattice/1000.xlsx'
#conf['datasets']['moments']['xlsx'][1001]='lattice/1001.xlsx'
#conf['datasets']['moments']['xlsx'][1002]='lattice/1002.xlsx'
#conf['datasets']['moments']['xlsx'][1002]='lattice/1002.xlsx'
#conf['datasets']['moments']['xlsx'][1003]='lattice/1003.xlsx'
#conf['datasets']['moments']['xlsx'][1004]='lattice/1004.xlsx'
#conf['datasets']['moments']['xlsx'][1005]='lattice/1005.xlsx'
#conf['datasets']['moments']['xlsx'][1006]='lattice/1006.xlsx'

#conf['datasets']['moments']['norm']={}
#for k in conf['datasets']['moments']['xlsx']: conf['datasets']['moments']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}

# AN

#conf['datasets']['AN']={}
#
#
#conf['datasets']['AN']['norm']={}
#
#conf['datasets']['AN']['filters']=[]
#
#
#conf['datasets']['AN']['xlsx']={}
#
#conf['datasets']['AN']['xlsx'][1000]='AN_pp/expdata/1000.xlsx' # BRAHMS pim 2.3
#conf['datasets']['AN']['xlsx'][1001]='AN_pp/expdata/1001.xlsx' # BRAHMS pim 4
#conf['datasets']['AN']['xlsx'][1002]='AN_pp/expdata/1002.xlsx' # BRAHMS pip 2.3
#conf['datasets']['AN']['xlsx'][1003]='AN_pp/expdata/1003.xlsx' # BRAHMS pip 4
#conf['datasets']['AN']['xlsx'][2000]='AN_pp/expdata/2000.xlsx' # STAR piz 04
#conf['datasets']['AN']['xlsx'][2001]='AN_pp/expdata/2001.xlsx' # STAR piz 3.3
#conf['datasets']['AN']['xlsx'][2002]='AN_pp/expdata/2002.xlsx' # STAR piz 3.68
#conf['datasets']['AN']['xlsx'][2003]='AN_pp/expdata/2003.xlsx' # STAR piz 3.7
#
#for k in conf['datasets']['AN']['xlsx']: conf['datasets']['AN']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}
