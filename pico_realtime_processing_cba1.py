"""
 pico_realtime_processing.py
 
  
 2017 CBA1 Mooring Realtime Data parsing and archiving.

 2017-07-16: copy original routine from 2017 SP03 and modify for 2017 deployment CBA1
 2017-03-31: use pandas for excel read instead of readXlsx()
 2016-12-12: add calculation to correct oxygen optode for salinity (aanderaa optodes
        have internal salinity set to 0 for basic operation... fresh water equivalent)
 2016-08-02: add class for met data from ketch - auto download and plot or choose from existing excel file
 2016-08-01: S.Bell
    cmocean colormap names updated (for V1.0)
"""

# Standard library.
import datetime
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
from scipy import interpolate
from netCDF4 import date2num, num2date

# Plotting Stack
import matplotlib as mpl
mpl.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, WeekdayLocator, MonthLocator, DayLocator, HourLocator, DateFormatter
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean

# User Stack
import calc.aanderaa_corrO2_sal as optode_O2_corr


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 05, 27)
__modified__ = datetime.datetime(2016, 05, 27)
__version__  = "0.1.1"
__status__   = "Development"

"""-------------------------- Initialization params -----------------------------------------"""


"""-------------------------------------- Main ----------------------------------------------"""

parser = argparse.ArgumentParser(description='RUDICS PICO Prawler Data File')
parser.add_argument('DataPath', metavar='DataPath', type=str,
               help='full path to file')


args = parser.parse_args()

### RUDICSead data as copied from the rudix server.
# http://ketch.pmel.noaa.gov/~taodata/view-picoprawl.html
data_dic = {}
startrow = ''
data_field = False
prw_cast_id = 1
with open(args.DataPath) as f:
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
            O2, O2_Temp, CO2_Temp, Phase, pCO2 = [], [], [], [], []
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
                                    'DO_Sat':DO_Sat,
                                    'O2':O2,
                                    'O2Temp':O2_Temp,
                                    'CO2Temp':CO2_Temp,
                                    'Phase':Phase,
                                    'pCO2':pCO2}]
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




