"""
Background:
===========
  prawler2netcdf_pico.py
 
 
Purpose:
========
  Script to convert cleaned-up radiometer data from RUDICS stream into netcdf file

  See extras/17ck-met_rad.py for cleanup script
 
 Compatibility:
 ==============
  python >=3.6 ? 
  python 2.7 - Tested

"""

# Standard library.
import datetime

# System Stack
import argparse

# Scientific stack.
import numpy as np
import pandas as pd
from netCDF4 import date2num, num2date

# User Stack
import io_utils.EcoFOCI_netCDF_write as EcF_write
import io_utils.ConfigParserLocal as ConfigParserLocal

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 12, 19)
__modified__ = datetime.datetime(2016, 12, 19)
__version__  = "0.1.1"
__status__   = "Development"


"""-------------------------------------- Main ----------------------------------------------"""

parser = argparse.ArgumentParser(description='RUDICS PICO Prawler Data File')
parser.add_argument('DataPath', metavar='DataPath', type=str,
               help='full path to file')
parser.add_argument('ConfigFile', metavar='ConfigFile', type=str,
               help='full path to nc config file')
parser.add_argument('OutPreFix', metavar='OutPreFix', type=str,
               help='prefix for output file')
parser.add_argument('-rads','--rads', action="store_true",
               help='radiometer data')


args = parser.parse_args()

### RUDICSead data as copied from the rudix server.
# http://ketch.pmel.noaa.gov/~taodata/view-picoprawl.html
data_dic = {}

df_rads = pd.read_csv(args.DataPath)
recnum = len(df_rads)

    
EPIC_VARS_dict = ConfigParserLocal.get_config(args.ConfigFile,'yaml')
data_dic['time'] = date2num([datetime.datetime.strptime(x,'%Y-%m-%d %H:%M:%S') for x in df_rads['timestamp'].values],'hours since 1900-01-01')
data_dic['SWRAD'] = df_rads[' rads_swr'].values
data_dic['LWRAD'] = df_rads[' rads_lwr'].values
data_dic['LWRAD_net'] = df_rads['net_lwr'].values
data_dic['CASETEMP'] = df_rads[' case temp (degC)'].values
data_dic['DOMETEMP'] = df_rads[' dome_temp (degC)'].values

if args.rads:
    #TODO: Need to make a profile_samplenum id for sequential ordering
    #create new netcdf file
    ncinstance = EcF_write.NetCDF_Create_Timeseries(savefile='data/' + args.OutPreFix  + '_rads' '.nc')
    ncinstance.file_create()
    ncinstance.sbeglobal_atts(raw_data_file=args.DataPath.split('/')[-1], 
        History='File Created.')
    ncinstance.dimension_init(time_len=recnum)
    ncinstance.variable_init(EPIC_VARS_dict,udunits_time_str='hours since 1900-01-01T00:00:00Z')
    ncinstance.add_coord_data(time=range(1,recnum+1))
    ncinstance.add_data(EPIC_VARS_dict,data_dic=data_dic,missing_values=np.nan)
    ncinstance.close()
