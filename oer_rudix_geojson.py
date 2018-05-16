#!/usr/bin/env python

"""
 Background:
 --------
 oer_rudix_geojson.py
 
 
 Purpose:
 --------
 dump rudix locations to geojson file

 History:
 --------


"""

#System Stack
import datetime
import argparse
from io import BytesIO
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

#science stack
import numpy as np

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 6, 01)
__modified__ = datetime.datetime(2016, 6, 01)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'CTD', 'SeaWater', 'Cruise', 'derivations','ketch','prawler','geojson'



"""----------------------- Class definitions --------------------------------------"""

class KetchMetData(object):
    r"""Download and parse pico met data from ketch.pmel.noaa.gov"""

    @staticmethod
    def get_data(time,station,server):
        r"""Download data from the PMEL engineering system.

        Parameters
        ----------
        time : datetime
            Start Date and time for which data should be downloaded

        Returns
        -------
        a file-like object from which to read the data
        """
        url = ('{server}/tao-bin/show_spurs2prawl?prawloption=gps&'
                'progid=pico&platid={station}&start={time:%Y}-{time:%m}-{time:%d}&end=&output=text&sensor=all').format(time=time,station=station,server=server)
        fobj = urlopen(url)
        data = fobj.read()

        # Since the archive text format is embedded in HTML, look for the <PRE> tags
        data_start = data.find(b'<PRE>')
        data_end = data.find(b'</PRE>', data_start)

        # Grab the stuff *between* the <PRE> tags -- 6 below is len('<PRE>\n')
        buf = data[data_start + 6:data_end]
        return BytesIO(buf.strip())

    @staticmethod
    def parse(fobj):
        r"""Parse Ketch Pico Met Data.

        This parses the particular tabular layout of met data for the PICO system.

        Parameters
        ----------
        fobj : file-like object
            The file-like object from which the data should be read. This needs to be set up
            to return bytes when read, not strings.

        Returns
        -------
        dict of information used 
        """

        # Parse the actual data, only grabbing the columns for date, time, lat, lon, quality
        cols = (0, 1, 2, 3, 4)

        # Skip 1 header lines -- 1 for '=====' 
        data = np.genfromtxt(fobj, names="GMTDate, GMTTime, Lat, Lon, quality", dtype=('S10', 'S8', 'f8', 'f8', 'f8'), usecols=cols, skip_header=2, unpack=True)
        kdatetime = [datetime.datetime.strptime(a+'  '+c,'%Y-%m-%d %H:%M:%S') for a,c in zip(data['GMTDate'],data['GMTTime'])]
        return dict(datetime=(kdatetime, 'datetime'), Lat=(data['Lat'], 'Deg'), Lon=(data['Lon'], 'Deg'),
                    quality=(data['quality'], ''))

"""--------------------------- Main ------------------------------------------------"""

parser = argparse.ArgumentParser(description='RUDICS PICO Prawler Location Data File')
parser.add_argument('-DataPath','--DataPath', type=str,
               help='full path to file')
parser.add_argument('-rudix','--rudix', action="store_true",
               help='flag to auto download rudix loc data')
parser.add_argument('-sid','--stationid', type=str, default='SP03',
               help='rudix station id (defaults to 2017 M2 deployment')
parser.add_argument('-sd','--start_date', type=str, default='2017-04-28',
               help='rudix station startdate (defaults to 2017 M2 deployment. Use yyyy-mm-dd')
parser.add_argument("-s",'--ServerName', type=str,
               default="http://yawl.pmel.noaa.gov",
               help='server name, eg. http://yawl.pmel.noaa.gov')
args = parser.parse_args()

###
#

if args.rudix:

    kid = KetchMetData.get_data(time=datetime.datetime.strptime(args.start_date,'%Y-%m-%d'),
        station=args.stationid,
        server=args.ServerName)
    data = KetchMetData.parse(kid)

    header = '{"type": "FeatureCollection","features": ['
    geojson = []
    for index, val in enumerate(data['Lat'][0]):
        geojson =geojson + ['{{"type": "Feature","geometry": {{"type": "Point", "coordinates":  [ {lon},{lat} ]}},"properties": {{"datetime":"{datetime}","quality":"{quality}"}}}}'.format(
                lat=data['Lat'][0][index],
                lon=data['Lon'][0][index],
                datetime=data['datetime'][0][index],
                quality=data['quality'][0][index])]

    geojson = ",".join(geojson)
    geojson = header + geojson + ']}'

    print geojson