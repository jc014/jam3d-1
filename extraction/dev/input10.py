conf={}

############################################################################
#mcsamp
conf['ncpu']=10
conf['nruns']=1
conf['factor']=3.0
conf['kappa']=1.5
conf['tol']=1e-10
conf['itmax']=int(1e7)
conf['block size']=10

############################################################################
# params

conf['params']={}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conf['params']['pdf']={}
conf['params']['pdf']['widths1_uv']  ={'value':    5.65394e-01,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths1_dv']  ={'value':    6.01982e-01,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths1_sea'] ={'value':    6.19623e-01,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths2_uv']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths2_dv']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['pdf']['widths2_sea'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conf['params']['ffpi']={}
conf['params']['ffpi']['widths1_fav']  ={'value':    1.16836e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths1_ufav'] ={'value':    1.36697e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths2_fav']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['ffpi']['widths2_ufav'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conf['params']['ffk']={}
conf['params']['ffk']['widths1_fav']  ={'value':    1.32455e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths1_ufav'] ={'value':    1.85629e-01,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_fav']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['ffk']['widths2_ufav'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conf['params']['transversity']={}
conf['params']['transversity']['widths1_uv']  ={'value':    7.16247e-01,'min': 0,'max':1,'fixed':False}
conf['params']['transversity']['widths1_dv']  ={'value':    7.16247e-01,'min': 0,'max':1,'fixed':'widths1_uv'}
conf['params']['transversity']['widths1_sea'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':False}
conf['params']['transversity']['widths2_uv']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['widths2_dv']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['widths2_sea'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

conf['params']['transversity']['u N0 1']  ={'value':    2.57906e-01,'min': 0,'max':10,'fixed':False}
conf['params']['transversity']['u a0 1']  ={'value':   -7.62366e-02,'min':-1,'max':5,'fixed':False}
conf['params']['transversity']['u b0 1']  ={'value':    1.40636e+00,'min':-1,'max':5,'fixed':False}
conf['params']['transversity']['u c0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['u d0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

conf['params']['transversity']['d N0 1']  ={'value':   -3.99037e-01,'min':-10,'max':1,'fixed':False}
conf['params']['transversity']['d a0 1']  ={'value':    6.27674e-02,'min':-1,'max':5,'fixed':False}
conf['params']['transversity']['d b0 1']  ={'value':    4.20148e-01,'min':-1,'max':5,'fixed':False}
conf['params']['transversity']['d c0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['d d0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

conf['params']['transversity']['s N0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['s a0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['s b0 1']  ={'value':    1.00000e+01,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['s c0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['transversity']['s d0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

conf['params']['transversity']['ub N0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s N0 1'}
conf['params']['transversity']['ub a0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s a0 1'}
conf['params']['transversity']['ub b0 1'] ={'value':    1.00000e+01,'min': 0,'max':1,'fixed':'s b0 1'}
conf['params']['transversity']['ub c0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s c0 1'}
conf['params']['transversity']['ub d0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s d0 1'}

conf['params']['transversity']['db N0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s N0 1'}
conf['params']['transversity']['db a0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s a0 1'}
conf['params']['transversity']['db b0 1'] ={'value':    1.00000e+01,'min': 0,'max':1,'fixed':'s b0 1'}
conf['params']['transversity']['db c0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s c0 1'}
conf['params']['transversity']['db d0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s d0 1'}

conf['params']['transversity']['sb N0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s N0 1'}
conf['params']['transversity']['sb a0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s a0 1'}
conf['params']['transversity']['sb b0 1'] ={'value':    1.00000e+01,'min': 0,'max':1,'fixed':'s b0 1'}
conf['params']['transversity']['sb c0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s c0 1'}
conf['params']['transversity']['sb d0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'s d0 1'}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conf['params']['collinspi']={}
conf['params']['collinspi']['widths1_fav']  ={'value':    3.78222e-02,'min': 0,'max':1,'fixed':False}
conf['params']['collinspi']['widths1_ufav'] ={'value':    1.23700e-01,'min': 0,'max':1,'fixed':False}
conf['params']['collinspi']['widths2_fav']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['collinspi']['widths2_ufav'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

# favored
conf['params']['collinspi']['u N0 1']  ={'value':    1.36506e+00,'min': 0,'max':10,'fixed':False}
conf['params']['collinspi']['u a0 1']  ={'value':   -5.62896e-01,'min':-1,'max':5,'fixed':False}
conf['params']['collinspi']['u b0 1']  ={'value':    2.31378e+00,'min':-1,'max':5,'fixed':False}
conf['params']['collinspi']['u c0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['collinspi']['u d0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

conf['params']['collinspi']['db N0 1'] ={'value':    1.36506e+00,'min': 0,'max':1,'fixed':'u N0 1'}
conf['params']['collinspi']['db a0 1'] ={'value':   -5.62896e-01,'min': 0,'max':1,'fixed':'u a0 1'}
conf['params']['collinspi']['db b0 1'] ={'value':    2.31378e+00,'min': 0,'max':1,'fixed':'u b0 1'}
conf['params']['collinspi']['db c0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'u c0 1'}
conf['params']['collinspi']['db d0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'u d0 1'}

# unfavored
conf['params']['collinspi']['d N0 1']  ={'value':   -1.76183e+00,'min':-10,'max':5,'fixed':False}
conf['params']['collinspi']['d a0 1']  ={'value':   -4.14970e-01,'min':-1,'max':5,'fixed':False}
conf['params']['collinspi']['d b0 1']  ={'value':    2.72143e+00,'min':-1,'max':5,'fixed':False}
conf['params']['collinspi']['d c0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
conf['params']['collinspi']['d d0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}

conf['params']['collinspi']['ub N0 1'] ={'value':   -1.76183e+00,'min': 0,'max':1,'fixed':'d N0 1'}
conf['params']['collinspi']['ub a0 1'] ={'value':   -4.14970e-01,'min': 0,'max':1,'fixed':'d a0 1'}
conf['params']['collinspi']['ub b0 1'] ={'value':    2.72143e+00,'min': 0,'max':1,'fixed':'d b0 1'}
conf['params']['collinspi']['ub c0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'d c0 1'}
conf['params']['collinspi']['ub d0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'d d0 1'}

conf['params']['collinspi']['s N0 1']  ={'value':   -1.76183e+00,'min': 0,'max':1,'fixed':'d N0 1'}
conf['params']['collinspi']['s a0 1']  ={'value':   -4.14970e-01,'min': 0,'max':1,'fixed':'d a0 1'}
conf['params']['collinspi']['s b0 1']  ={'value':    2.72143e+00,'min': 0,'max':1,'fixed':'d b0 1'}
conf['params']['collinspi']['s c0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'d c0 1'}
conf['params']['collinspi']['s d0 1']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'d d0 1'}

conf['params']['collinspi']['sb N0 1'] ={'value':   -1.76183e+00,'min': 0,'max':1,'fixed':'d N0 1'}
conf['params']['collinspi']['sb a0 1'] ={'value':   -4.14970e-01,'min': 0,'max':1,'fixed':'d a0 1'}
conf['params']['collinspi']['sb b0 1'] ={'value':    2.72143e+00,'min': 0,'max':1,'fixed':'d b0 1'}
conf['params']['collinspi']['sb c0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'d c0 1'}
conf['params']['collinspi']['sb d0 1'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':'d d0 1'}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#conf['params']['collinsk']={}
#conf['params']['collinsk']['widths1_fav']  ={'value':   -1.25354e-03,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['widths1_ufav'] ={'value':    4.26794e-07,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['widths2_fav']  ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
#conf['params']['collinsk']['widths2_ufav'] ={'value':    0.00000e+00,'min': 0,'max':1,'fixed':True}
#
#conf['params']['collinsk']['u N0 1']  ={'value':    1,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['u a0 1']  ={'value': -0.5,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['u b0 1']  ={'value':    3,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['u c0 1']  ={'value':    0,'min': 0,'max':1,'fixed':True}
#conf['params']['collinsk']['u d0 1']  ={'value':    0,'min': 0,'max':1,'fixed':True}
#
#conf['params']['collinsk']['sb N0 1'] ={'value':    1,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['sb a0 1'] ={'value': -0.5,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['sb b0 1'] ={'value':    3,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['sb c0 1'] ={'value':    0,'min': 0,'max':1,'fixed':True}
#conf['params']['collinsk']['sb d0 1'] ={'value':    0,'min': 0,'max':1,'fixed':True}
#
#conf['params']['collinsk']['d N0 1']  ={'value':    1,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['d a0 1']  ={'value': -0.5,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['d b0 1']  ={'value':    3,'min': 0,'max':1,'fixed':False}
#conf['params']['collinsk']['d c0 1']  ={'value':    0,'min': 0,'max':1,'fixed':True}
#conf['params']['collinsk']['d d0 1']  ={'value':    0,'min': 0,'max':1,'fixed':True}
#
#conf['params']['collinsk']['db N0 1'] ={'value':    1,'min': 0,'max':1,'fixed':'d N0 1'}
#conf['params']['collinsk']['db a0 1'] ={'value': -0.5,'min': 0,'max':1,'fixed':'d a0 1'}
#conf['params']['collinsk']['db b0 1'] ={'value':    3,'min': 0,'max':1,'fixed':'d b0 1'}
#conf['params']['collinsk']['db c0 1'] ={'value':    0,'min': 0,'max':1,'fixed':'d c0 1'}
#conf['params']['collinsk']['db d0 1'] ={'value':    0,'min': 0,'max':1,'fixed':'d d0 1'}
#
#conf['params']['collinsk']['ub N0 1'] ={'value':    1,'min': 0,'max':1,'fixed':'d N0 1'}
#conf['params']['collinsk']['ub a0 1'] ={'value': -0.5,'min': 0,'max':1,'fixed':'d a0 1'}
#conf['params']['collinsk']['ub b0 1'] ={'value':    3,'min': 0,'max':1,'fixed':'d b0 1'}
#conf['params']['collinsk']['ub c0 1'] ={'value':    0,'min': 0,'max':1,'fixed':'d c0 1'}
#conf['params']['collinsk']['ub d0 1'] ={'value':    0,'min': 0,'max':1,'fixed':'d d0 1'}
#
#conf['params']['collinsk']['s N0 1']  ={'value':    1,'min': 0,'max':1,'fixed':'d N0 1'}
#conf['params']['collinsk']['s a0 1']  ={'value': -0.5,'min': 0,'max':1,'fixed':'d a0 1'}
#conf['params']['collinsk']['s b0 1']  ={'value':    3,'min': 0,'max':1,'fixed':'d b0 1'}
#conf['params']['collinsk']['s c0 1']  ={'value':    0,'min': 0,'max':1,'fixed':'d c0 1'}
#conf['params']['collinsk']['s d0 1']  ={'value':    0,'min': 0,'max':1,'fixed':'d d0 1'}



############################################################################
# set data sets

conf['datasets']={}
conf['datasets']['sidis']={}
conf['datasets']['sidis']['xlsx']={}

# upol
#conf['datasets']['sidis']['xlsx'][1000]='sidis/expdata/1000.xlsx'  # |  proton   | pi+   | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1001]='sidis/expdata/1001.xlsx'  # |  proton   | pi-   | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1004]='sidis/expdata/1004.xlsx'  # |  deuteron | pi+   | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx'  # |  deuteron | pi-   | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1002]='sidis/expdata/1002.xlsx'  # |  proton   | k+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1003]='sidis/expdata/1003.xlsx'  # |  proton   | k-    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx'  # |  deuteron | k+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1007]='sidis/expdata/1007.xlsx'  # |  deuteron | k-    | M_Hermes | hermes 

# transversity
conf['datasets']['sidis']['xlsx'][3000]='sidis/expdata/3000.xlsx' # | proton   | pi+    | AUTcollins       | hermes     | x                    
conf['datasets']['sidis']['xlsx'][3003]='sidis/expdata/3003.xlsx' # | proton   | pi+    | AUTcollins       | hermes     | z                    
conf['datasets']['sidis']['xlsx'][3026]='sidis/expdata/3026.xlsx' # | proton   | pi+    | AUTcollins       | hermes     | pt                   
conf['datasets']['sidis']['xlsx'][3004]='sidis/expdata/3004.xlsx' # | proton   | pi-    | AUTcollins       | hermes     | x                    
conf['datasets']['sidis']['xlsx'][3018]='sidis/expdata/3018.xlsx' # | proton   | pi-    | AUTcollins       | hermes     | z                    
conf['datasets']['sidis']['xlsx'][3016]='sidis/expdata/3016.xlsx' # | proton   | pi-    | AUTcollins       | hermes     | pt                   
#conf['datasets']['sidis']['xlsx'][3006]='sidis/expdata/3006.xlsx' # | proton   | pi0    | AUTcollins       | hermes     | z                    
#conf['datasets']['sidis']['xlsx'][3014]='sidis/expdata/3014.xlsx' # | proton   | pi0    | AUTcollins       | hermes     | x                    
#conf['datasets']['sidis']['xlsx'][3015]='sidis/expdata/3015.xlsx' # | proton   | pi0    | AUTcollins       | hermes     | pt                   

conf['datasets']['sidis']['xlsx'][3025]='sidis/expdata/3025.xlsx' # | proton   | pi+    | AUTcollins       | compass    | x                    
conf['datasets']['sidis']['xlsx'][3010]='sidis/expdata/3010.xlsx' # | proton   | pi+    | AUTcollins       | compass    | z                    
conf['datasets']['sidis']['xlsx'][3027]='sidis/expdata/3027.xlsx' # | proton   | pi+    | AUTcollins       | compass    | pt                   
conf['datasets']['sidis']['xlsx'][3005]='sidis/expdata/3005.xlsx' # | proton   | pi-    | AUTcollins       | compass    | x                    
conf['datasets']['sidis']['xlsx'][3013]='sidis/expdata/3013.xlsx' # | proton   | pi-    | AUTcollins       | compass    | z                    
conf['datasets']['sidis']['xlsx'][3012]='sidis/expdata/3012.xlsx' # | proton   | pi-    | AUTcollins       | compass    | pt                   

conf['datasets']['sidis']['xlsx'][4000]='sidis/expdata/4000.xlsx' # | deuteron | pi+    | AUTcollins       | compass    | x                    
conf['datasets']['sidis']['xlsx'][4002]='sidis/expdata/4002.xlsx' # | deuteron | pi+    | AUTcollins       | compass    | z                    
conf['datasets']['sidis']['xlsx'][4001]='sidis/expdata/4001.xlsx' # | deuteron | pi+    | AUTcollins       | compass    | pt                   
conf['datasets']['sidis']['xlsx'][4003]='sidis/expdata/4003.xlsx' # | deuteron | pi-    | AUTcollins       | compass    | x                    
conf['datasets']['sidis']['xlsx'][4005]='sidis/expdata/4005.xlsx' # | deuteron | pi-    | AUTcollins       | compass    | z                    
conf['datasets']['sidis']['xlsx'][4004]='sidis/expdata/4004.xlsx' # | deuteron | pi-    | AUTcollins       | compass    | pt                   

#conf['datasets']['sidis']['xlsx'][3001]='sidis/expdata/3001.xlsx' # | neutron  | pi-    | AUTcollins       | jlab       | x                    
#conf['datasets']['sidis']['xlsx'][3002]='sidis/expdata/3002.xlsx' # | neutron  | pi+    | AUTcollins       | jlab       | x                  

#conf['datasets']['sidis']['xlsx'][3007]='sidis/expdata/3007.xlsx' # | proton   | k+     | AUTcollins       | hermes     | x                    
#conf['datasets']['sidis']['xlsx'][3008]='sidis/expdata/3008.xlsx' # | proton   | k+     | AUTcollins       | hermes     | z                    
#conf['datasets']['sidis']['xlsx'][3024]='sidis/expdata/3024.xlsx' # | proton   | k+     | AUTcollins       | hermes     | pt                   
#conf['datasets']['sidis']['xlsx'][3017]='sidis/expdata/3017.xlsx' # | proton   | k-     | AUTcollins       | hermes     | x                    
#conf['datasets']['sidis']['xlsx'][3023]='sidis/expdata/3023.xlsx' # | proton   | k-     | AUTcollins       | hermes     | z                    
#conf['datasets']['sidis']['xlsx'][3021]='sidis/expdata/3021.xlsx' # | proton   | k-     | AUTcollins       | hermes     | pt                   
#
#
#conf['datasets']['sidis']['xlsx'][6000]='sidis/expdata/6000.xlsx' # | proton   | k-     | AUTcollins       | compass    | pt                   
#conf['datasets']['sidis']['xlsx'][6001]='sidis/expdata/6001.xlsx' # | proton   | k-     | AUTcollins       | compass    | x                    
#conf['datasets']['sidis']['xlsx'][6002]='sidis/expdata/6002.xlsx' # | proton   | k-     | AUTcollins       | compass    | z                    
#conf['datasets']['sidis']['xlsx'][6003]='sidis/expdata/6003.xlsx' # | proton   | k+     | AUTcollins       | compass    | pt                   
#conf['datasets']['sidis']['xlsx'][6004]='sidis/expdata/6004.xlsx' # | proton   | k+     | AUTcollins       | compass    | x                    
#conf['datasets']['sidis']['xlsx'][6005]='sidis/expdata/6005.xlsx' # | proton   | k+     | AUTcollins       | compass    | z                    
#
#
#conf['datasets']['sidis']['xlsx'][4006]='sidis/expdata/4006.xlsx' # | deuteron | k+     | AUTcollins       | compass    | x                    
#conf['datasets']['sidis']['xlsx'][4008]='sidis/expdata/4008.xlsx' # | deuteron | k+     | AUTcollins       | compass    | z                    
#conf['datasets']['sidis']['xlsx'][4007]='sidis/expdata/4007.xlsx' # | deuteron | k+     | AUTcollins       | compass    | pt                   
#conf['datasets']['sidis']['xlsx'][4009]='sidis/expdata/4009.xlsx' # | deuteron | k-     | AUTcollins       | compass    | x                    
#conf['datasets']['sidis']['xlsx'][4011]='sidis/expdata/4011.xlsx' # | deuteron | k-     | AUTcollins       | compass    | z                    
#conf['datasets']['sidis']['xlsx'][4010]='sidis/expdata/4010.xlsx' # | deuteron | k-     | AUTcollins       | compass    | pt                   



conf['datasets']['sidis']['norm']={}
for idx in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][idx]={'value':1,'fixed':True,'min':0,'max':1} 

conf['datasets']['sidis']['filters']={}
for idx in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['filters'][idx]="z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9"







