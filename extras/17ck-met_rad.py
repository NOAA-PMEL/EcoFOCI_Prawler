#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 08:07:40 2018

@author: bell
"""
import datetime


#%% Load initial raw data and format output
# Step 1

fp = '/Volumes/WDC_internal/Users/bell/ecoraid/2017/Moorings/17ckitaepr2a/raw/met/chukchi_met.txt'

param='LWR'
with open(fp) as f:
    cont=False
    for line in f:
        if (param in line):
            lineparts=line.split()
            try:
                dt=datetime.datetime.strptime(lineparts[1]+' '+lineparts[2],'%m/%d/%Y %H:%M:%S')
                cont=True
                count=0
            except ValueError:
                continue
        if '\n' == line:
            cont=False
        if not (param in line) and cont :
           t1 = ",".join(line.split()[0:5])
           t2 = ",".join(line.split()[5:])
           print(dt+datetime.timedelta(minutes=count),',',t1.replace('\n',''))
           count +=1
           print(dt+datetime.timedelta(minutes=count),',',t2.replace('\n',''))
           count +=1


#%% Load formatted output
# Step 2           

import pandas as pd

fp_lwr = '/Volumes/WDC_internal/Users/bell/ecoraid/2017/Moorings/17ckitaepr2a/rawconverted/met/17ckitaepr2a_lwr.txt'
fp_swr = '/Volumes/WDC_internal/Users/bell/ecoraid/2017/Moorings/17ckitaepr2a/rawconverted/met/17ckitaepr2a_swr.txt'

df_lwr = pd.read_csv(fp_lwr,parse_dates=True,index_col=['timestamp'])
df_swr = pd.read_csv(fp_swr,parse_dates=True,index_col=['timestamp'])

df_rads = df_swr.join(df_lwr, how='outer', lsuffix='_swr', rsuffix='_lwr')

df_rads_1min = df_rads.resample('1t').median()

#%% Calculate LWR Net

sigma = 5.6704e-8
df_rads_1min['net_lwr'] = (sigma*((273.15+df_rads_1min[' dome_temp (degC)'])**4)) + df_rads_1min[' rads_lwr'] \
                            - 3.5*sigma*(((273.15+df_rads_1min[' dome_temp (degC)'])**4)-((273.15+df_rads_1min[' case temp (degC)'])**4))
                            
df_rads_1min.to_csv('/Volumes/WDC_internal/Users/bell/ecoraid/2017/Moorings/17ckitaepr2a/rawconverted/met/17ckitaepr2a_rads.txt')