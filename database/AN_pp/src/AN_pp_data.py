# -*- coding: utf-8 -*-
"""
AN_pp_data.py - program to convert AN pp .dat data files to .xlsx
"""

import os
import numpy as np
import pandas as pd

# Import data as Excel files

files = os.listdir('./')
files = [f for f in files if f.endswith('.dat')]

names = [f.strip('.dat') for f in files]

names_dict = {"pim.2.3.200": '1000', "pim.4.200": '1001', "pip.2.3.200": '1002', "pip.4.200": '1003',
              "piz.04": '2000', "piz.3.3": '2001', "piz.3.68": '2002', "piz.3.7": '2003'}

had_dict = {"pim.2.3.200": 'pi-', "pim.4.200": 'pi-', "pip.2.3.200": 'pi+', "pip.4.200": 'pi+',
            "piz.04": 'pi0', "piz.3.3": 'pi0', "piz.3.68": 'pi0', "piz.3.7": 'pi0'}

for i in range(len(files)):
    df = pd.read_table(files[i], sep='\s+', names=['xF',
                                                   'value', 'stat_err_u', 'sys_err_u', 'pT'])
    had_col = [had_dict[names[i]] for j in range(len(df.index))]
    b_col = ['BRAHMS' for j in range(len(df.index))]
    s_col = ['STAR' for j in range(len(df.index))]
    obs_col = ['AN' for j in range(len(df.index))]
    tar_col = ['p' for j in range(len(df.index))]
    rs_col = [200.0 for j in range(len(df.index))]
    df['rs'] = rs_col
    df['target'] = tar_col
    df['hadron'] = had_col
    if 'pi+' in had_col:
        df['col'] = b_col
    elif 'pi-' in had_col:
        df['col'] = b_col
    elif 'pi0' in had_col:
        df['col'] = s_col
    df['obs'] = obs_col
    writer = pd.ExcelWriter(names_dict[names[i]] + r'.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, sheet_name='Sheet1')
    writer.save()
