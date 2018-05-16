"""
 prawler2netcdf_pico.py
 
 Description:
    Using text output from Ketch data server (which converts the rudix data), create
        an archival format for ERDDAP and other dissemination.

 2016 ITAE Mooring Realtime Data parsing and archiveing (from pico).

 2018-05-09: add ability to save wx data into netcdf from pico format
 2016-12-12: add calculation to correct oxygen optode for salinity (aanderaa optodes
        have internal salinity set to 0 for basic operation... fresh water equivalent)
"""

# Standard library.
import datetime

# System Stack
import argparse

# Scientific stack.
import numpy as np
import seawater as sw
from netCDF4 import date2num, num2date

# User Stack
import calc.aanderaa_corrO2_sal as optode_O2_corr
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
parser.add_argument('-met','--met', action="store_true",
               help='meteorological data')
parser.add_argument('-is2D','--is2D', action="store_true",
               help='1D ragged arrays')
parser.add_argument('-perdive','--perdive', action="store_true",
               help='file per downcast')
args = parser.parse_args()

### RUDICSead data as copied from the rudix server.
# http://ketch.pmel.noaa.gov/~taodata/view-picoprawl.html
data_dic = {}
startrow = ''
data_field = False
prw_cast_id = 1
recnum = 0
with open(args.DataPath) as f:
    if not args.met:
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
                sample, time, Depth, Temp, Cond, Salinity, DO, DO_Temp, Chl, Turb, SigmaT, DO_Sat = [], [], [], [], [], [], [], [], [], [], [], []
        

            if (len(line) == 0) and (startrow != ''):
                startrow = ''
                data_dic[prw_cast_id] = [{'sample':sample, 
                                        'time':time,
                                        'Depth':Depth, 
                                        'Temperature':Temp, 
                                        'Cond':Cond, 
                                        'Salinity':Salinity, 
                                        'Oxy_Conc':DO, 
                                        'Oxy_Temperature':DO_Temp, 
                                        'Chlorophyll':Chl, 
                                        'Turbidity':Turb,
                                        'SigmaT':SigmaT,
                                        'Oxy_Sat':DO_Sat}]
                prw_cast_id += 1            

            if (k >= startrow) and (startrow != '') and data_field:
                sample = sample + [line.strip().split()[0]]
                time = time + [date2num(datetime.datetime.strptime(line.strip().split()[1]+' '+line.strip().split()[2],
                        '%Y-%m-%d %H:%M:%S'),
                        'hours since 1900-01-01T00:00:00Z')]
                Depth = Depth + [np.float(line.strip().split()[3])]
                Temp = Temp + [np.float(line.strip().split()[4])]
                Cond = Cond + [np.float(line.strip().split()[5])]
                Salinity = Salinity + [np.float(line.strip().split()[6])]
                #DO = DO + [np.float(line.strip().split()[7])]
                DO_Temp = DO_Temp + [np.float(line.strip().split()[8])]
                Chl = Chl + [np.float(line.strip().split()[9])]
                Turb = Turb + [np.float(line.strip().split()[10])]
                # calculate sigmaT at 0db gauge pressure (s, t, p=0)
                SigmaT = SigmaT + [sw.eos80.dens0(s=np.float(line.strip().split()[6]),t=np.float(line.strip().split()[4]))-1000.]

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
                recnum += 1
                
                # remove first entry for files from html/wget routines with 
                try:
                    data_dic.pop(1)
                except:
                    pass
                try:
                    data_dic.pop(2)
                except:
                    pass

    elif args.met:
        for k, line in enumerate(f.readlines()):
            line = line.strip()

            if 'Suspect' in line:
                continue

            if '<PRE>' in line:
                data_field = True
            if '</PRE>' in line:
                data_field = False
                line = []

            if '=========' in line and data_field:  # get start line of data
                startrow = k + 1
                sample, time, uwind, vwind, wspd, wdir, rh, at, bp = [], [], [], [], [], [], [], [], []
        

            if (len(line) == 0) and (startrow != ''):
                startrow = ''
                data_dic['met'] = {
                                        'time':time,
                                        'eastward_wind':uwind, 
                                        'northward_wind':vwind, 
                                        'wind_speed':wspd, 
                                        'wind_from_direction':wdir, 
                                        'relative_humidity':rh, 
                                        'air_temperature':at, 
                                        'air_pressure_at_sealevel':bp}
                prw_cast_id += 1            

            if (k >= startrow) and (startrow != '') and data_field:
                time = time + [date2num(datetime.datetime.strptime(line.strip().split()[0]+' '+line.strip().split()[1],
                        '%Y-%m-%d %H:%M:%S'),
                        'hours since 1900-01-01T00:00:00Z')]
                uwind = uwind + [np.float(line.strip().split()[2])]
                vwind = vwind + [np.float(line.strip().split()[3])]
                wspd = wspd + [np.float(line.strip().split()[4])]
                wdir = wdir + [np.float(line.strip().split()[5])]
                rh = rh + [np.float(line.strip().split()[6])]
                at = at + [np.float(line.strip().split()[7])]
                bp = bp + [np.float(line.strip().split()[8])]
                recnum += 1


    
EPIC_VARS_dict = ConfigParserLocal.get_config(args.ConfigFile,'yaml')

if args.met:
    #TODO: Need to make a profile_samplenum id for sequential ordering
    #create new netcdf file
    ncinstance = EcF_write.NetCDF_Create_Timeseries(savefile='data/' + args.OutPreFix  + '_met' '.nc')
    ncinstance.file_create()
    ncinstance.sbeglobal_atts(raw_data_file=args.DataPath.split('/')[-1], 
        History='File Created.')
    ncinstance.dimension_init(time_len=recnum)
    ncinstance.variable_init(EPIC_VARS_dict,udunits_time_str='hours since 1900-1-1')
    ncinstance.add_coord_data(time=range(1,recnum+1))
    ncinstance.add_data(EPIC_VARS_dict,data_dic=data_dic['met'],missing_values=np.nan)
    ncinstance.close()

if args.is2D:

    #find max obs in a profile to set for all samples
    obs = []
    for key in data_dic.keys():
        obs = obs + [len(data_dic[key][0]['sample'])]
    #create new netcdf file
    ncinstance = EcF_write.NetCDF_Create_Profile_Ragged2D(savefile='data/' + args.OutPreFix + '.nc')
    ncinstance.file_create()
    ncinstance.sbeglobal_atts(raw_data_file=args.DataPath.split('/')[-1], 
        History='File Created.  Aanderaa Optode Dissolved O2 compensated for Salinity/Depth')
    ncinstance.dimension_init(profilenum_len=len(data_dic.keys()),obsnum_len=np.max(obs))
    ncinstance.variable_init(EPIC_VARS_dict,udunits_time_str='hours since 1900-1-1')
    ncinstance.add_coord_data(profile_num=data_dic.keys(), obs_num=range(1,np.max(obs)+1))
    for ind, profile_num, in enumerate(data_dic.keys()):
        print "Adding Profile {profile_num}".format(profile_num=profile_num)
        ncinstance.add_data(EPIC_VARS_dict,
                            profile_num=ind, 
                            data_dic=data_dic[profile_num][0],
                            missing_values=np.nan)
    ncinstance.close()

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
            print "Adding Profile {profile_num}".format(profile_num=key)
            ncinstance.add_data(EPIC_VARS_dict,
                                data_dic=data_dic[key][0],
                                missing_values=np.nan)

            ncinstance.close()
        except KeyError:
            pass