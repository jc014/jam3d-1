conf={}

#--setups
conf['bootstrap']=False
conf['flat par']=True
conf['ftol']=1e-8
conf['ncpus']=4
conf['Q20']   = 1.27**2
conf['order'] = 'LO'

#--datasets

conf['datasets']={}

conf['datasets']['sidis']={}
conf['datasets']['sidis']['filters']=["z>0.2 and z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9"]
conf['datasets']['sidis']['xlsx']={}

#--Unpol multiplicities (HERMES only)
conf['datasets']['sidis']['xlsx'][1000]='sidis/expdata/1000.xlsx'  # |  proton   | pi+   | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1001]='sidis/expdata/1001.xlsx'  # |  proton   | pi-   | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1004]='sidis/expdata/1004.xlsx'  # |  deuteron | pi+   | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx'  # |  deuteron | pi-   | M_Hermes | hermes

#--Collins effect SIDIS
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
conf["datasets"]["sidis"]["xlsx"][3027]="sidis/expdata/3027.xlsx"  #  compass    proton   pi+ pt
conf["datasets"]["sidis"]["xlsx"][3025]="sidis/expdata/3025.xlsx"  #  compass    proton  pi+   x
conf["datasets"]["sidis"]["xlsx"][3010]="sidis/expdata/3010.xlsx"  #  compass    proton  pi+  z
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

conf['datasets']['sidis']['norm']={}
for k in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][k]={'value':1.00000e+00,'fixed':False,'min':8.00000e-01,'max':1.20000e+00}


#--parameters
conf['params']={}


#--Parameters in gaussian approximation, parton model:
#--TMD PDF:
conf['params']['pdf']={}
conf['params']['pdf']['widths1_uv']  ={'value':    5.40344e-01,'min': 0.1,'max':0.8,'fixed':False}
conf['params']['pdf']['widths2_uv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['pdf']['widths1_dv']  ={'value':    5.40344e-01,'min': 0.1,'max':0.8,'fixed':'widths1_uv'}
conf['params']['pdf']['widths2_dv']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}
conf['params']['pdf']['widths1_sea'] ={'value':    5.85373e-01,'min': 0.1,'max':0.8,'fixed':False}
conf['params']['pdf']['widths2_sea'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_uv'}

#--TMD FF:
conf['params']['ffpi']={}
conf['params']['ffpi']['widths1_fav']  ={'value':    1.22011e-01,'min': 0.05,'max':0.3,'fixed':False}
conf['params']['ffpi']['widths2_fav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffpi']['widths1_ufav'] ={'value':    1.43665e-01,'min': 0.05,'max':0.3,'fixed':False}
conf['params']['ffpi']['widths2_ufav'] ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}

conf['params']['ffk']={}
conf['params']['ffk']['widths1_fav']   ={'value':    1.32333e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffk']['widths1_ufav']  ={'value':    2.02660e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}

conf['params']['ffh']={}
conf['params']['ffh']['widths1_fav']   ={'value':    1.24050e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_fav']   ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':True}
conf['params']['ffh']['widths1_ufav']  ={'value':    1.43730e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffh']['widths2_ufav']  ={'value':    0.00000e+00,'min':-1,'max':1,'fixed':'widths2_fav'}


# Transversity
conf['params']['transversity+']={}
conf['params']['transversity+']['widths1_uv']  ={'value':    7.11165e-01, 'min':0.2, 'max':0.8, 'fixed': False}
conf['params']['transversity+']['widths1_dv']  ={'value':    2.12697e-01, 'min':0.2, 'max':0.8, 'fixed': 'widths1_uv'}
conf['params']['transversity+']['widths1_sea'] ={'value':    5.35004e-01, 'min':0, 'max':1, 'fixed': True}
conf['params']['transversity+']['widths2_uv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': True}
conf['params']['transversity+']['widths2_dv']  ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}
conf['params']['transversity+']['widths2_sea'] ={'value':    0.00000e+00, 'min':0, 'max':1, 'fixed': 'widths2_uv'}

conf['params']['transversity+']['g1 N']    ={'value':    3.87592e-01, 'min':  None, 'max':  None, 'fixed': False }
conf['params']['transversity+']['g1 a']    ={'value':   -6.23068169e-01, 'min':  -1.9, 'max':     1, 'fixed': False}
conf['params']['transversity+']['g1 b']    ={'value':    9.25741583e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['uv1 N']   ={'value':    3.47549e-01, 'min':  None, 'max':  None, 'fixed': False }
conf['params']['transversity+']['uv1 a']   ={'value':   -1.21835956e-01, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['transversity+']['uv1 b']   ={'value':    3.20766744e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['dv1 N']   ={'value':    1.52089e-01, 'min':  None, 'max':  None, 'fixed': False }
conf['params']['transversity+']['dv1 a']   ={'value':   -2.39874967e-01, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['transversity+']['dv1 b']   ={'value':    3.83902620e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['db1 N']   ={'value':    3.67609928e-02, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['transversity+']['db1 a']   ={'value':   -8.41360631e-01, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['transversity+']['db1 b']   ={'value':    5.31285539e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['ub1 N']   ={'value':    1.95464789e-02, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['transversity+']['ub1 a']   ={'value':    -9.93659187e-01, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['transversity+']['ub1 b']   ={'value':    8.38905814e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['s1 N']    ={'value':    0.00000e+00, 'min':     0, 'max':     1, 'fixed': False }
conf['params']['transversity+']['s1 a']    ={'value':   1.34706224e-01, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['transversity+']['s1 b']    ={'value':    6.00759596e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['sb1 N']   ={'value':    7.46109845e-07, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['transversity+']['sb1 a']   ={'value':   3.83495317e-01, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['transversity+']['sb1 b']   ={'value':    4.61209808e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['sea1 N']  ={'value':    5.71081196e-03, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['transversity+']['sea1 a']  ={'value':   -1.36329697e+00, 'min':  -1.9, 'max':    -1, 'fixed': False}
conf['params']['transversity+']['sea1 b']  ={'value':    4.74721050e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['transversity+']['sea2 b']  ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': 'sea1 b'}




# Collins TMD FF:
conf['params']['collinspi+']={}

conf['params']['collinspi+']['widths1_fav']   ={'value':    1.22288e-01,'min':0.005, 'max':0.5, 'fixed':False}
conf['params']['collinspi+']['widths1_ufav'] ={'value':    6.11536e-02,'min':0.005, 'max':0.5, 'fixed':False}
conf['params']['collinspi+']['widths2_fav']   ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':True}
conf['params']['collinspi+']['widths2_ufav'] ={'value':    0.00000e+00,'min':0, 'max':1, 'fixed':'widths2_fav'}

conf['params']['collinspi+']['g1 N']  ={'value':    2.95437e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['collinspi+']['g1 a']  ={'value':    1.00469e+00, 'min': -1.8 , 'max':     2,'fixed':False}
conf['params']['collinspi+']['g1 b']  ={'value':    6.85766e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                
conf['params']['collinspi+']['u1 N']  ={'value':    2.67821e-02, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['collinspi+']['u1 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':     2,'fixed':False}
conf['params']['collinspi+']['u1 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['collinspi+']['d1 N']  ={'value':    2.99974e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['collinspi+']['d1 a']  ={'value':   -6.89477e-01, 'min': -1.8 , 'max':     2,'fixed':False}
conf['params']['collinspi+']['d1 b']  ={'value':    4.79992e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['collinspi+']['s1 N']  ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['collinspi+']['s1 a']  ={'value':    3.00305e-01, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
conf['params']['collinspi+']['s1 b']  ={'value':    1.83178e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                   
conf['params']['collinspi+']['c1 N']  ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['collinspi+']['c1 a']  ={'value':   -5.05798e-02, 'min': -1.8 , 'max':     2,'fixed':False}
conf['params']['collinspi+']['c1 b']  ={'value':    3.19952e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['collinspi+']['b1 N']  ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['collinspi+']['b1 a']  ={'value':   -1.59541e+00, 'min': -1.8 , 'max':     2,'fixed':False}
conf['params']['collinspi+']['b1 b']  ={'value':    4.50102e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['collinspi+']['ub1 N'] ={'value':    2.99974e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['collinspi+']['ub1 a'] ={'value':   -6.89477e-01, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
conf['params']['collinspi+']['ub1 b'] ={'value':    4.79992e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                  
conf['params']['collinspi+']['db1 N'] ={'value':    2.67821e-02, 'min':  0.0 , 'max':     1,'fixed':'u1 N'}
conf['params']['collinspi+']['db1 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':     2,'fixed':'u1 a'}
conf['params']['collinspi+']['db1 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    10,'fixed':'u1 b'}
                                                                                                  
conf['params']['collinspi+']['sb1 N'] ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['collinspi+']['sb1 a'] ={'value':    3.00305e-01, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
conf['params']['collinspi+']['sb1 b'] ={'value':    1.83178e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                  
conf['params']['collinspi+']['cb1 N'] ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':'c1 N'}
conf['params']['collinspi+']['cb1 a'] ={'value':   -5.05798e-02, 'min': -1.8 , 'max':     2,'fixed':'c1 a'}
conf['params']['collinspi+']['cb1 b'] ={'value':    3.19952e+00, 'min':  0   , 'max':    10,'fixed':'c1 b'}
                                                                                                  
conf['params']['collinspi+']['bb1 N'] ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':'b1 N'}
conf['params']['collinspi+']['bb1 a'] ={'value':   -1.59541e+00, 'min': -1.8 , 'max':     2,'fixed':'b1 a'}
conf['params']['collinspi+']['bb1 b'] ={'value':    4.50102e+00, 'min':  0   , 'max':    10,'fixed':'b1 b'}




