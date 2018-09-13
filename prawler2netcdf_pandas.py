
"""
 Background:
 ===========
 prawler2netcdf_pandas.py
 
 
 Purpose:
 ========
 Modification of prawler2netcdf_pico.py to make same formatted output using pandas dataframe


 History:
 ========

 Compatibility:
 ==============
 python >=3.6 
 python 2.7 

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
__created__  = datetime.datetime(2018, 12, 19)
__modified__ = datetime.datetime(2018, 12, 19)
__version__  = "0.1.0"
__status__   = "Development"


"""-------------------------------------- Main ----------------------------------------------"""

parser = argparse.ArgumentParser(description='Pandas PICO Prawler Data File')
parser.add_argument('DataPath', metavar='DataPath', type=str,
               help='full path to file')
parser.add_argument('ConfigFile', metavar='ConfigFile', type=str,
               help='full path to nc config file')
parser.add_argument('OutPreFix', metavar='OutPreFix', type=str,
               help='prefix for output file')
parser.add_argument('-perdive','--perdive', action="store_true",
               help='file per downcast')
args = parser.parse_args()

df = pd.read_csv(args.DataPath,parse_dates=True)
df['time_num'] = [date2num(datetime.datetime.strptime(x,'%Y-%m-%d %H:%M:%S'),
                        'hours since 1900-01-01T00:00:00Z') for x in df.datetime]    
EPIC_VARS_dict = ConfigParserLocal.get_config(args.ConfigFile,'yaml')

data_dic = {}

dfg = df.groupby('divenumber')
for prw_cast_id in dfg.groups.keys():
    sample, time, Depth, Temp, Cond, Salinity = [], [], [], [], [], [] 
    DO, DO_Temp, Chl, Turb, SigmaT, DO_Sat = [], [], [], [], [], []

    data_dic[prw_cast_id] = [{'sample':dfg.get_group(prw_cast_id).index.values, 
                            'time':dfg.get_group(prw_cast_id).time_num.values,
                            'Depth':dfg.get_group(prw_cast_id).press.values, 
                            'Temperature':dfg.get_group(prw_cast_id).temp.values, 
                            'Cond':dfg.get_group(prw_cast_id).cond.values, 
                            'Salinity':dfg.get_group(prw_cast_id).sal.values, 
                            'Oxy_Conc':dfg.get_group(prw_cast_id).DO.values, 
                            'Oxy_Temperature':dfg.get_group(prw_cast_id).c2.values, 
                            'Chlorophyll':dfg.get_group(prw_cast_id).Chlor.values, 
                            'Turbidity':dfg.get_group(prw_cast_id).Turb.values,
                            'SigmaT':dfg.get_group(prw_cast_id).SigmaT.values,
                            'Oxy_Sat':dfg.get_group(prw_cast_id).DO_Sat.values}]


if args.perdive:
    #find max obs in a profile to set for all samples
    for key in data_dic.keys():
        try:
            obs = len(data_dic[key][0]['sample'])
            #create new netcdf file
            ncinstance = EcF_write.NetCDF_Create_Profile_Ragged1D(savefile='data/' +args.OutPreFix  + '_p' + str(key).zfill(4) + '.nc')
            ncinstance.file_create()
            ncinstance.sbeglobal_atts(raw_data_file=args.DataPath.split('/')[-1], 
                History='File Created.  Aanderaa Optode Dissolved O2 compensated for Salinity/Depth')
            ncinstance.dimension_init(recnum_len=obs)
            ncinstance.variable_init(EPIC_VARS_dict)
            ncinstance.add_coord_data(range(1,np.max(obs)+1))
            print("Adding Profile {profile_num}".format(profile_num=key))
            ncinstance.add_data(EPIC_VARS_dict,
                                data_dic=data_dic[key][0],
                                missing_values=np.nan)

            ncinstance.close()
        except KeyError:
            pass