"""
 pico_profile_processing.py
 
  
 2017 SP03 Mooring Realtime Data parsing and archiving.

 2017-08-08: Port cross section routine to make profile plots

"""

# Standard library.
import datetime
import os
from io import BytesIO
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
    
# System Stack
import argparse

# Scientific stack.
import numpy as np
import pandas as pd
import seawater as sw

# Plotting Stack
import matplotlib as mpl
mpl.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cmocean

# User Stack
import calc.aanderaa_corrO2_sal as optode_O2_corr
from plots.profile_plot import CTDProfilePlot


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 05, 27)
__modified__ = datetime.datetime(2016, 05, 27)
__version__  = "0.1.1"
__status__   = "Development"

"""-------------------------- Quick Subroutines -----------------------------------------"""


def castdirection(depth):
    """determin index of upcast and downcast"""
    downcast = [0,np.argmax(depth)+1]
    upcast = [np.argmax(depth)+1,len(depth)]

    return (downcast,upcast)

"""-------------------------------------- Main ----------------------------------------------"""
parser = argparse.ArgumentParser(description='RUDICS PICO Prawler Data File')
parser.add_argument('DataPath', metavar='DataPath', type=str,
               help='full path to file')
parser.add_argument('latlon', metavar='DataPath', type=str, nargs=2,
               help='lat lon')

args = parser.parse_args()



### RUDICSead data as copied from the rudix server.
# http://ketch.pmel.noaa.gov/~taodata/view-picoprawl.html
data_dic = {}
startrow = ''
data_field = False
prw_cast_id = 1
with open(args.DataPath) as f:
    ### TODO: Replace with class that either ffetches from ketch or a datafile
    for k, line in enumerate(f.readlines()):
        line = line.strip()

        if 'Suspect' in line:
            continue

        if '<PRE>' in line:
            data_field = True
        if '</PRE>' in line:
            data_field = False
            line = []

        if '>' in line and data_field:  # get start line of data
            startrow = k + 1
            sample, Date, Time, Depth, Temp, Cond, Salinity, DO, DO_Temp, Chl, Turb, SigmaT, DO_Sat = [], [], [], [], [], [], [], [], [], [], [], [], []

        if (len(line) == 0) and (startrow != ''):
            startrow = ''
            data_dic[prw_cast_id] = [{'sample':sample, 
                                    'Date':Date, 'Time':Time,
                                    'Depth':Depth, 
                                    'Temp':Temp, 
                                    'Cond':Cond, 
                                    'Salinity':Salinity, 
                                    'DO':DO, 
                                    'DO_Temp':DO_Temp, 
                                    'Chl':Chl, 
                                    'Turb':Turb,
                                    'SigmaT':SigmaT,
                                    'DO_Sat':DO_Sat}]
        prw_cast_id += 1            

        if (k >= startrow) and (startrow != '') and data_field:
            sample = sample + [line.strip().split()[0]]
            Date = Date + [line.strip().split()[1]]
            Time = Time + [line.strip().split()[2]]
            Depth = Depth + [np.float(line.strip().split()[3])]
            Temp = Temp + [np.float(line.strip().split()[4])]
            Cond = Cond + [np.float(line.strip().split()[5])]
            Salinity = Salinity + [np.float(line.strip().split()[6])]
            #DO = DO + [np.float(line.strip().split()[7])]
            DO_Temp = DO_Temp + [np.float(line.strip().split()[8])]
            Chl = Chl + [np.float(line.strip().split()[9])]
            Turb = Turb + [np.float(line.strip().split()[10])]
            # calculate sigmaT at 0db gauge pressure (s, t, p=0)
            SigmaT = SigmaT + [sw.eos80.dens0(np.float(line.strip().split()[6]),np.float(line.strip().split()[4]))-1000.]

            # apply salinity and depth corrections to oxygen optode and recalc percentsat
            O2_corr = optode_O2_corr.O2_dep_comp(oxygen_conc=np.float(line.strip().split()[7]),
                                                 depth=np.float(line.strip().split()[3]))
            O2_corr = optode_O2_corr.O2_sal_comp(oxygen_conc=O2_corr,
                                                 salinity=np.float(line.strip().split()[6]),
                                                 temperature=np.float(line.strip().split()[4]))
            DO = DO + [optode_O2_corr.O2_molar2umkg(oxygen_conc=O2_corr,
                                             salinity=np.float(line.strip().split()[6]),
                                             temperature=np.float(line.strip().split()[4]),
                                             pressure=np.float(line.strip().split()[3]))]
            DO_Sat = DO_Sat + [optode_O2_corr.O2PercentSat(oxygen_conc=O2_corr, 
                                             temperature=np.float(line.strip().split()[4]),
                                             salinity=np.float(line.strip().split()[6]),
                                             pressure=np.float(line.strip().split()[3]))]

# remove first entry for files from html/wget routines with 
#data_dic.pop(1)
#data_dic.pop(2)

count = 1
for k in sorted(data_dic.keys()):

    #######################
    #
    # Plots
    if not data_dic[k][0]['Depth']:
        print("Cast {profile} is empty".format(profile=count))
        continue
    else:
        print("Cast {profile} is being processed".format(profile=count))

    PrawlerPlot = CTDProfilePlot()
    downInd,upInd = castdirection(np.array(data_dic[k][0]['Depth']))

    depth = np.array(data_dic[k][0]['Depth'])
    Temp = np.array(data_dic[k][0]['Temp'])
    DO_Temp = np.array(data_dic[k][0]['DO_Temp'])
    Salinity = np.array(data_dic[k][0]['Salinity'])

    if not os.path.exists('images/profile_' + str(count) ):
        os.makedirs('images/profile_' + str(count))

    ########## CTD
    ### temperature
    (plt, fig) = PrawlerPlot.plot1plot(epic_key=['T_28','T2_35'],
                     xdata=[Temp,DO_Temp],
                     ydata=[depth],
                     xlabel='Temperature (C)',
                     updown=['d','d'],
                     maxdepth=np.max(depth),
                     plotpoints=True)

    ptitle = PrawlerPlot.add_title(cruiseid='',
                      fileid=args.DataPath.split('/')[-1],
                      castid=str(count),
                      castdate=datetime.datetime.strptime(data_dic[k][0]['Date'][0] + ' ' + data_dic[k][0]['Time'][0],'%Y-%m-%d %H:%M:%S'),
                      lat=args.latlon[0],
                      lon=args.latlon[1])

    t = fig.suptitle(ptitle)
    t.set_y(1.06)
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0], DefaultSize[1]*2) )

    plt.savefig('images/profile_'+str(count)+'/profile_'+str(count)+'_temperature.png', bbox_inches='tight', dpi = (300))
    plt.close()

    ### salinity
    (plt, fig) = PrawlerPlot.plot1plot(epic_key=['S_41'],
                     xdata=[Salinity],
                     ydata=[depth],
                     xlabel='Temperature (C)',
                     updown=['d'],
                     maxdepth=np.max(depth),
                     plotpoints=True)

    ptitle = PrawlerPlot.add_title(cruiseid='',
                      fileid=args.DataPath.split('/')[-1],
                      castid=str(count),
                      castdate=datetime.datetime.strptime(data_dic[k][0]['Date'][0] + ' ' + data_dic[k][0]['Time'][0],'%Y-%m-%d %H:%M:%S'),
                      lat=args.latlon[0],
                      lon=args.latlon[1])

    t = fig.suptitle(ptitle)
    t.set_y(1.06)
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0], DefaultSize[1]*2) )

    plt.savefig('images/profile_'+str(count)+'/profile_'+str(count)+'_salinity.png', bbox_inches='tight', dpi = (300))
    plt.close()

    ### temp,sal
    (plt, fig) = PrawlerPlot.plot2var(epic_key=['T_28','','S_41'],
                     xdata=[Temp,None,Salinity],
                     ydata=depth,
                     xlabel=['Temperature (C)','','Salinity (PSU)'],
                     maxdepth=np.max(depth))

    ptitle = PrawlerPlot.add_title(cruiseid='',
                      fileid=args.DataPath.split('/')[-1],
                      castid=str(count),
                      castdate=datetime.datetime.strptime(data_dic[k][0]['Date'][0] + ' ' + data_dic[k][0]['Time'][0],'%Y-%m-%d %H:%M:%S'),
                      lat=args.latlon[0],
                      lon=args.latlon[1])

    t = fig.suptitle(ptitle)
    t.set_y(1.06)
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0], DefaultSize[1]*2) )

    plt.savefig('images/profile_'+str(count)+'/profile_'+str(count)+'_TS.png', bbox_inches='tight', dpi = (300))
    plt.close()

    count+=1

