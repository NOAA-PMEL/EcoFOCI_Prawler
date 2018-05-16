#!/usr/bin/env python

"""
oer_rudix_processing.py
 
 Purpose:
 --------

 Ingest and organize Rudics transmitted data

 Discussion:  
 -----------
 Found on systems such as the Prawlers (profiling crawlers)

 Each instrument is passed as a chunk of data with a tag and list of lines to expect.
 Data needs to have calibrations or characterizations applied

 Example data with CTD/Aanderaa/Wetlabs

CTD 09/23/2017 18:01:56  08 046
3.600000 7.380500 3.223540 

AADI 09/23/2017 18:01:56  08 046
350.942 7.498 31.740 31.740 39.518 7.778 822.000 772.300 613.000

WETL 09/23/2017 18:01:56  08 046
0695 0183 0700 0053 0566

 History:
 --------
 2018-03-27: Begin wholesale modification of routines.  This script began as an interim
    way to evaluate data coming from Scott Stalin.


"""

from __future__ import absolute_import

#System Stack
import warnings
import datetime
import os
import sys
import argparse

# Scientific stack.
import numpy as np
import pandas as pd

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 06, 07)
__modified__ = datetime.datetime(2018, 06, 07)
__version__  = "0.2.0"
__status__   = "Development"
__keywords__ = 'along track','csv','timeseries','wavegliders','prawler'

__all__ = ['RUDIX_SERVICE_Profile',]

"""------------------------------- Data Read ----------------------------------------"""



class RUDIX_SERVICE_Profile(pd.DataFrame):
    r"""
    Translate profile data coming from the rudics service (usually provided by PMEL Engineering)
     and is often the raw way Prawler data is transmitted.

    Each data burst is a complete profile from each of the instruments of interest
    """
    def __init__(self, data=None, index=None, columns=None, 
                 divenum=0, longitude=None, latitude=None, missing=1e35,
                 dtype=None, copy=False):

      super(RUDIX_SERVICE_Profile, self).__init__(data=data, index=index,
                                columns=columns, dtype=dtype,
                                copy=copy)
      self.divenum = divenum
      self.missing = missing
      self.latitude = latitude
      self.longitude = longitude

    def __reduce__(self):
        return self.__class__, (
            pd.DataFrame(self), # NOTE Using that type(data)==DataFrame and the
                              # the rest of the arguments of DataFrame.__init__
                              # to defaults, the constructors acts as a
                              # copy constructor.
            None,
            None,
            self.divenum,
            self.missing,
            self.longitude,
            self.latitude,
            None,
            False,
        )


def RUDIX_CTD(fid,nrows=0):
    r"""
    assuming data with the following characteristics:

    pressure temperature conductivity
    (dbar) (degreeC) (Sm/?)

    """
    columns=['pressure','temperature','conductivity']
    return(pd.read_table(fid, header=None, index_col=None, names=columns, nrows=nrows, delim_whitespace=True))


def RUDIX_AADI(fid,nrows=0):
    r"""
    assuming data with the following characteristics:

    c1 c2 c3 c4 c5 c6 c7 c8 c9

    """
    columns=['c1','c2','c3','c4','c5','c6','c7','c8','c9']
    return(pd.read_table(fid, header=None, index_col=None, names=columns, nrows=nrows, delim_whitespace=True))


def RUDIX_WETL(fid,nrows):
    r"""
    assuming data with the following characteristics:

    channel_1 counts_1 channel_2 counts_2 counts_3

    channels are the wevelength of the transmission, the counts_3 channel is 
        often associated with temperature.

    """
    columns=['channel_1','counts_1','channel_2','counts_2','counts_3']
    return(pd.read_table(fid, header=None, index_col=None, names=columns, nrows=nrows, delim_whitespace=True))


"""------------------------------- lat/lon ----------------------------------------"""

def latlon_convert(Mooring_Lat, Mooring_Lon):
    
    tlat = Mooring_Lat.strip().split() #deg min dir
    lat = float(tlat[0]) + float(tlat[1]) / 60.
    if tlat[2] == 'S':
        lat = -1 * lat
        
    tlon = Mooring_Lon.strip().split() #deg min dir
    lon = float(tlon[0]) + float(tlon[1]) / 60.
    if tlon[2] == 'E':
        lon = -1 * lon
        
    return (lat, lon)

"""------------------------------- MAIN ----------------------------------------"""

parser = argparse.ArgumentParser(description='Rudix Data Processing')
parser.add_argument('DataPath', metavar='DataPath', type=str, help='full path to file')
parser.add_argument('-nc', '--netcdf_out',
    action="store_true", 
    help='save each profile as netcdf')

args = parser.parse_args()

profile_count = 0

with open(args.DataPath) as f:

    for k, line in enumerate(f.readlines()):
        try:
            line_array = line.strip().split()
        except:
            continue
        try:
            line_array[0]
        except:
            continue

        if line_array[0] == 'CTD': #start of a CTD record
            time_start = datetime.datetime.strptime(line_array[1]+' '+line_array[2],'%m/%d/%Y %H:%M:%S')
            no_samples = int(line_array[4])
            CTD = RUDIX_CTD(f,no_samples)
            print(CTD)

        elif line_array[0] == 'AADI': #end CTD record
            time_start = datetime.datetime.strptime(line_array[1]+' '+line_array[2],'%m/%d/%Y %H:%M:%S')
            no_samples = int(line_array[4])
            AADI = RUDIX_AADI(f,no_samples)
            print(AADI)

        elif line_array[0] == 'WETL': #start of wetlabs record
            time_start = datetime.datetime.strptime(line_array[1]+' '+line_array[2],'%m/%d/%Y %H:%M:%S')
            no_samples = int(line_array[4])
            WETL = RUDIX_WETL(f,no_samples)
            print(WETL)

        else:
            continue






