"""
prawler_wget.py

Purpose:
	Connect to yawl.pmel.noaa.gov (engineering run system) to retrieve prawler data from 2016 ITAE Bering Sea Mooring

 2018-05-02: S.Bell - with ketch being replaced by yawl (servers http for https)
 	make server name an argument
 2017-04-26: Make routines more generic by adding a few more arguments
"""
#System Stack
import datetime
import argparse

import wget

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 6, 01)
__modified__ = datetime.datetime(2016, 6, 01)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'CTD', 'SeaWater', 'Cruise', 'derivations','yawl','prawler'

"""--------------------------------helper Routines---------------------------------------"""


"""--------------------------------main Routines---------------------------------------"""

parser = argparse.ArgumentParser(description='yawl Data Retrieval')
parser.add_argument('Project', metavar='Project', type=str,
               help='project name e.g. 2016_ITAE')
parser.add_argument('PlatformID', metavar='PlatformID', type=str,
               help='yawl platform ID e.g. PITW')
parser.add_argument("-s",'--ServerName', type=str,
			   default="http://yawl.pmel.noaa.gov",
               help='server name, eg. http://yawl.pmel.noaa.gov')
parser.add_argument("-a",'--all', action="store_true", help='grabs all 2016 prawler data and places it in to one file')
parser.add_argument("-as",'--all_since', type=str, help='grabs all 2016 prawler data since specified (yyyy-mm-dd) and places it in to one file')
parser.add_argument("-ad",'--all_daily', action="store_true", help='grabs all 2016 prawler data and places it in to daily files')
parser.add_argument("-t",'--today', action="store_true", help='grabs todays prawler day and places it in a file')
parser.add_argument("-ud",'--userday', type=str, help='choose day to retrieve -- format yyyy-mm-dd')
parser.add_argument("-met",'--met', action="store_true", 
	help ='retrieve all met data')
parser.add_argument("-loc",'--loc', action="store_true", 
	help='retrieve all loc data')
          
args = parser.parse_args()

# retrieve all data - one file
if args.all:
	filename = args.Project+'_all.download'
	url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=ctd&progid=PICO&platid='+args.PlatformID+ \
		    '&start=2017-04-25&end=&selectall=true&output=text&zgrid=0&background=white&colormap=ocean_r&ncol=6'
	wget.download(url, filename, bar=wget.bar_thermometer)

# retrieve all data - one file
if args.met:
	filename = args.Project+'_all.met.download'
	url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=met&progid=pico&platid='+args.PlatformID+ \
		    '&start=2017-04-25&end=&selectall=true&output=text&zgrid=0&background=white&colormap=ocean_r&ncol=6'
	wget.download(url, filename, bar=wget.bar_thermometer)

if args.loc:
	filename = args.Project+'_all.loc.download'
	url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=gps&progid=pico&platid='+args.PlatformID+ \
		    '&start=2017-04-01&end=&output=text'
	wget.download(url, filename, bar=wget.bar_thermometer)


# retrieve all data - multiple files
if args.all_daily:
	startday = datetime.datetime(2016, 05, 04)
	today = datetime.datetime.now()
	dates = [startday + datetime.timedelta(x) for x in range(0,365,1)]
	for i,day in enumerate(dates):
		if not day > today:
			print "Retrieving {0}".format(day.strftime('%Y%m%d'))
			filename = args.Project + '_daily.'+ day.strftime('%Y%m%d')+'.download'
			url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=ctd&progid=PICO&platid='+args.PlatformID+ \
				    '&start=' + day.strftime('%Y-%m-%d') + '&end=' + (day + datetime.timedelta(1)).strftime('%Y-%m-%d') + \
				    '&selectall=true&output=text&zgrid=0&background=white&colormap=ocean_r&ncol=6'
			wget.download(url, filename, bar=False)
		else:
			break

# retrieve one day of data
if args.today:
	today = datetime.datetime.now()
	tomorrow = today + datetime.timedelta(1)
	filename = args.Project+'_daily.'+today.strftime('%Y%m%d')+'.download'
	url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=ctd&progid=PICO&platid='+args.PlatformID+ \
		    '&start=' + today.strftime('%Y-%m-%d') + '&end=' + tomorrow.strftime('%Y-%m-%d') + \
		    '&selectall=true&output=text&zgrid=0&background=white&colormap=ocean_r&ncol=6'
	wget.download(url, filename, bar=False)

# retrieve one day of data - users choice
if args.userday:
	today = datetime.datetime.strptime(args.userday, "%Y-%m-%d")
	tomorrow = today + datetime.timedelta(1)
	filename = args.Project+'_daily.'+today.strftime('%Y%m%d')+'.download'
	url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=ctd&progid=PICO&platid='+args.PlatformID+ \
		    '&start=' + today.strftime('%Y-%m-%d') + '&end=' + tomorrow.strftime('%Y-%m-%d') + \
		    '&selectall=true&output=text&zgrid=0&background=white&colormap=ocean_r&ncol=6'
	wget.download(url, filename, bar=False)

# retrieve all data - one file , specifying the start date
if args.all_since:
	filename = args.Project+'.pt2.download'
	url = args.ServerName+'/tao-bin/show_spurs2prawl?prawloption=ctd&progid=PICO&platid='+args.PlatformID+ \
		    '&start='+args.all_since+'&end=&selectall=true&output=text&zgrid=0&background=white&colormap=ocean_r&ncol=6'
	wget.download(url, filename, bar=wget.bar_thermometer)

