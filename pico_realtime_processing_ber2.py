"""
 pico_realtime_processing.py
 
  
 2017 SP03 Mooring Realtime Data parsing and archiving.

 2018-05-01: copy original routine from 2017 SP03 for 2018 BER2 deployment
    ketch was replaced by yawl (which supports https) 
2017-04-26: copy original routine from 2016 ITAE and modify for 2017 deployment SP03
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

### some mpl specif settings for fonts and plot style
mpl.rcParams['svg.fonttype'] = 'none'
### Plot Data
mpl.rcParams['font.size'] = 6.
mpl.rcParams['axes.labelsize'] = 8.
mpl.rcParams['xtick.labelsize'] = 10.
mpl.rcParams['ytick.labelsize'] = 8.

### Plot Data
mpl.rcParams['font.size'] = 6.
mpl.rcParams['axes.labelsize'] = 8.
mpl.rcParams['xtick.labelsize'] = 10.
mpl.rcParams['ytick.labelsize'] = 8.

"""-------------------------------------- subroutines ---------------------------------------"""

def cmap_map(function,cmap):
    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(map( reduced_cmap, step_list))
    new_LUT = np.array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return mpl.colors.LinearSegmentedColormap('colormap',cdict,1024)

def RHtoTd( RH, Tc):
    """ 
    Parameters
    ----------
    RH : array_like
         Relative Humidity (expressed in percent)
    Tc : array-like
         Temperature (in Degrees Celsius)
    
    Returns
    -------
    Outputs : array_like
              Dewpoint Temperature (in Degrees Celsius)
        
    Notes
    -----
        Uses:     Es=6.11*10.0**(7.5*Tc/(237.7+Tc))
                  E=(RH*Es)/100
                  Tdc=(-430.22+237.7*np.log(E))/(-1.0*np.log(E)+19.08)
    
    Examples
    --------
    TODO
    
    References
    ----------
    TODO
    
    Modifications
    -------------
    """

    RH = np.array(RH)
    Tc = np.array(Tc)
    
    Es=6.11*10.0**(7.5*Tc/(237.7+Tc))

    E=(RH*Es)/100

    Tdc=(-430.22+237.7*np.log(E))/(-1.0*np.log(E)+19.08)

    return(Tdc)


"""----------------------- Class definitions --------------------------------------"""

class KetchMetData(object):
    r"""Download and parse pico met data from ketch.pmel.noaa.gov"""

    @staticmethod
    def get_data(time):
        r"""Download data from the PMEL engineering system.

        Parameters
        ----------
        time : datetime
            Start Date and time for which data should be downloaded

        Returns
        -------
        a file-like object from which to read the data
        """
        url = ('https://yawl.pmel.noaa.gov/tao-bin/show_spurs2prawl?prawloption=met&progid=pico&'
                'platid=BER2&start={time:%Y}-{time:%m}-{time:%d}&end=&output=text&sensor=all').format(time=time)
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

        # Parse the actual data, only grabbing the columns for date, time, RH, AT, and wind speed
        cols = (0, 1, 4, 6, 7)

        # Skip 1 header lines -- 1 for '=====' 
        data = np.genfromtxt(fobj, names="GMTDate, GMTTime, Spd, RH, AT", dtype=('S10', 'S8', 'f8', 'f8', 'f8'), usecols=cols, skip_header=3, unpack=True)
        kdatetime = [datetime.datetime.strptime(a+'  '+c,'%Y-%m-%d %H:%M:%S') for a,c in zip(data['GMTDate'],data['GMTTime'])]
        return dict(datetime=(kdatetime, 'datetime'), AT=(data['AT'], 'DegC'), RH=(data['RH'], '%'),
                    wind=(data['Spd'], 'm/s'))





"""-------------------------------------- Main ----------------------------------------------"""

parser = argparse.ArgumentParser(description='RUDICS PICO Prawler Data File')
parser.add_argument('DataPath', metavar='DataPath', type=str,
               help='full path to file')
parser.add_argument('-met_rudix','--met_rudix', action="store_true",
               help='flag to auto download rudix met data')
parser.add_argument('-met','--WindTempDataPath', type=str,
               help='optional full path to met filefile')
parser.add_argument('-fg','--FillGaps', action="store_true",
               help='Interpolate and Fill Gaps in bin averaged data')
parser.add_argument('-grid','--Gridded', action="store_true",
               help='put on regular grid')
parser.add_argument('-i','--image', action="store_true",
               help='Make an image plot')
parser.add_argument('-csv_out','--csv_out', type=str,
               help='save chosen parameter to file')


args = parser.parse_args()

tw_time,tw_temp,tw_wind,tw_td = {},{},{},{}
if args.WindTempDataPath:
    tempwind_file = args.WindTempDataPath
    print "Reading file {0}".format(tempwind_file)
    W = pd.read_excel(tempwind_file)
    for index, row in W.iterrows():
        tw_time[index] = pd.to_datetime(row['Date Time']).to_pydatetime()
    tw_time = tw_time.values()
    tw_temp = W['AT ( C)'].values
    tw_td   = RHtoTd(W['RH'],W['AT ( C)'])
    tw_wind = W['Wind Spd (m/s)'].values

if args.met_rudix:
    print "\nDownloading Met data from Rudix system"
    kid = KetchMetData.get_data(datetime.datetime(2018,4,30))
    data = KetchMetData.parse(kid)
    for index, val in enumerate(data['AT'][0]):
        tw_time[index] = data['datetime'][0][index]
        tw_temp[index] = val
        tw_td[index] = RHtoTd(data['RH'][0][index],data['AT'][0][index])
        tw_wind[index] = data['wind'][0][index]
    tw_time = tw_time.values()
    tw_temp = tw_temp.values()
    tw_td = tw_td.values()
    tw_wind = tw_wind.values()


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
try:
    data_dic.pop(min(data_dic.keys()))
except:
    pass

### vertically grid data to evenly space gridspoints
interval = 0.25
press_grid = np.arange(0,50,interval) #1m

mesh_grid_s, mesh_grid_t, mesh_grid_o, mesh_grid_chl = [], [], [], []
mesh_grid_sig, mesh_grid_turb, mesh_grid_osat, mesh_grid_stats = [], [], [], []
date_time = []

for k in data_dic.keys():

    mesh_depth_s, mesh_depth_t, mesh_depth_o, mesh_depth_chl = [], [], [], []
    mesh_depth_sig, mesh_depth_turb, mesh_depth_osat, mesh_depth_stats = [], [], [], []
    irreg_depth = np.array(data_dic[k][0]['Depth'])
    irreg_sal   = np.array(data_dic[k][0]['Salinity'])
    irreg_temp  = np.array(data_dic[k][0]['Temp'])
    irreg_oxy   = np.array(data_dic[k][0]['DO'])
    irreg_osat  = np.array(data_dic[k][0]['DO_Sat'])
    irreg_chlor = np.array(data_dic[k][0]['Chl'])
    irreg_sigmat= np.array(data_dic[k][0]['SigmaT'])
    irreg_turb  = np.array(data_dic[k][0]['Turb'])
    cast_date   = mpl.dates.date2num(datetime.datetime.strptime(data_dic[k][0]['Date'][0] + ' ' + data_dic[k][0]['Time'][0],'%Y-%m-%d %H:%M:%S'))
    
    #TODO: update with groupby statement
    for pg in press_grid:
        ireg_ind = np.where((irreg_depth > pg) & (irreg_depth <= pg+interval))
        mesh_depth_s = np.hstack((mesh_depth_s, np.median(irreg_sal[ireg_ind])))
        mesh_depth_sig = np.hstack((mesh_depth_sig, np.median(irreg_sigmat[ireg_ind])))
        mesh_depth_t = np.hstack((mesh_depth_t, np.median(irreg_temp[ireg_ind])))
        mesh_depth_o = np.hstack((mesh_depth_o, np.median(irreg_oxy[ireg_ind])))
        mesh_depth_osat = np.hstack((mesh_depth_osat, np.median(irreg_osat[ireg_ind])))
        mesh_depth_chl = np.hstack((mesh_depth_chl, np.median(irreg_chlor[ireg_ind])))
        mesh_depth_turb = np.hstack((mesh_depth_turb, np.median(irreg_turb[ireg_ind])))
        mesh_depth_stats = np.hstack((mesh_depth_stats, ireg_ind[0].size))
    
    if args.FillGaps:
        mask = np.isnan(mesh_depth_s)
        mesh_depth_s[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_s[~mask], right=-100000)
        mask = np.isnan(mesh_depth_t)
        mesh_depth_t[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_t[~mask], right=-100000)
        mask = np.isnan(mesh_depth_o)
        try:
            mesh_depth_o[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_o[~mask], right=-100000)
        except ValueError: #handles samples with all nan's
            mesh_depth_o[0]  = 0.0
            mesh_depth_o[-1] = 0.0
            mask = np.isnan(mesh_depth_o)
            mesh_depth_o[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_o[~mask], right=-100000)
        mask = np.isnan(mesh_depth_osat)
        try:
            mesh_depth_osat[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_osat[~mask], right=-100000)
        except ValueError: #handles samples with all nan's
            mesh_depth_osat[0]  = 0.0
            mesh_depth_osat[-1] = 0.0
            mask = np.isnan(mesh_depth_osat)
            mesh_depth_osat[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_osat[~mask], right=-100000)
        mask = np.isnan(mesh_depth_chl)
        try:
            mesh_depth_chl[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_chl[~mask], right=-100000)
        except ValueError: #handles samples with all nan's
            mesh_depth_chl[0]  = 0.0
            mesh_depth_chl[-1] = 0.0
            mask = np.isnan(mesh_depth_chl)
            mesh_depth_chl[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_chl[~mask], right=-100000)
        mask = np.isnan(mesh_depth_turb)
        try:
            mesh_depth_turb[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_turb[~mask], right=-100000)
        except ValueError: #handles samples with all nan's
            mesh_depth_turb[0]  = 0.0
            mesh_depth_turb[-1] = 0.0
            mask = np.isnan(mesh_depth_turb)
            mesh_depth_turb[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_turb[~mask], right=-1000)
        mask = np.isnan(mesh_depth_sig)
        mesh_depth_sig[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_sig[~mask], right=-1000)        
        mesh_depth_stats[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), mesh_depth_stats[~mask], right=-1000)        

    date_time = date_time + [cast_date]

    mesh_grid_s = mesh_grid_s + [mesh_depth_s]
    mesh_grid_t = mesh_grid_t + [mesh_depth_t]
    mesh_grid_o = mesh_grid_o + [mesh_depth_o]
    mesh_grid_osat = mesh_grid_osat + [mesh_depth_osat]
    mesh_grid_chl = mesh_grid_chl + [mesh_depth_chl]
    mesh_grid_sig = mesh_grid_sig + [mesh_depth_sig]
    mesh_grid_turb = mesh_grid_turb + [mesh_depth_turb]
    mesh_grid_stats = mesh_grid_stats + [mesh_depth_stats]

date_timetemp=date_time
date_time = np.array(date_time)

#grid time data
if args.Gridded:
    # put data on a regular grid (needed for imageshow() but not for contourf )
    dt = 1.0/24.0
    tmp = num2date(date_time.min(),'Days since 0001-1-1')
    dt_min = date2num(tmp -datetime.timedelta(seconds=60*tmp.minute + tmp.second),'Days since 0001-1-1')
    time_grid = np.arange(dt_min,date_time.max(),dt) #grid limits -> set to top of hour
    #grid_bounds = np.meshgrid(time_grid,press_grid)
    mesh_grid_sf = interpolate.interp2d(press_grid,date_time,mesh_grid_s)
    mesh_grid_s = mesh_grid_sf(press_grid,time_grid)
    mesh_grid_tf = interpolate.interp2d(press_grid,date_time,mesh_grid_t)
    mesh_grid_t = mesh_grid_tf(press_grid,time_grid)
    mesh_grid_of = interpolate.interp2d(press_grid,date_time,mesh_grid_o)
    mesh_grid_o = mesh_grid_of(press_grid,time_grid)
    mesh_grid_osatf = interpolate.interp2d(press_grid,date_time,mesh_grid_osat)
    mesh_grid_osat = mesh_grid_osatf(press_grid,time_grid)
    mesh_grid_chlf = interpolate.interp2d(press_grid,date_time,mesh_grid_chl)
    mesh_grid_chl = mesh_grid_chlf(press_grid,time_grid)
    mesh_grid_sigf = interpolate.interp2d(press_grid,date_time,mesh_grid_sig)
    mesh_grid_sig = mesh_grid_sigf(press_grid,time_grid)
    mesh_grid_turbf = interpolate.interp2d(press_grid,date_time,mesh_grid_turb)
    mesh_grid_turb = mesh_grid_turbf(press_grid,time_grid)
    mesh_grid_statsf = interpolate.interp2d(press_grid,date_time,mesh_grid_stats)
    mesh_grid_stats = mesh_grid_statsf(press_grid,time_grid)

    date_time = time_grid
    
    #fill known bad data points with missing values
    #fill known bad data points with missing values
    ### may3 - missing instrument
    bad_times = []
    
    for bad_time in bad_times:
        time_range =  np.where((time_grid >= bad_time[0]) & (time_grid <= bad_time[1]))
        mesh_grid_s[time_range,:] = mesh_grid_s[time_range,:]*np.nan
        mesh_grid_t[time_range,:] = mesh_grid_t[time_range,:]*np.nan
        mesh_grid_o[time_range,:] = mesh_grid_o[time_range,:]*np.nan
        mesh_grid_osat[time_range,:] = mesh_grid_osat[time_range,:]*np.nan
        mesh_grid_chl[time_range,:] = mesh_grid_chl[time_range,:]*np.nan
        mesh_grid_sig[time_range,:] = mesh_grid_sig[time_range,:]*np.nan
        mesh_grid_turb[time_range,:] = mesh_grid_turb[time_range,:]*np.nan
        mesh_grid_stats[time_range,:] = mesh_grid_stats[time_range,:]*np.nan  

    ##top 3m
    bad_depth =  np.where(press_grid <= 3)
    mesh_grid_s[:,bad_depth] = mesh_grid_s[:,bad_depth]*np.nan
    mesh_grid_t[:,bad_depth] = mesh_grid_t[:,bad_depth]*np.nan
    mesh_grid_o[:,bad_depth] = mesh_grid_o[:,bad_depth]*np.nan
    mesh_grid_osat[:,bad_depth] = mesh_grid_osat[:,bad_depth]*np.nan
    mesh_grid_chl[:,bad_depth] = mesh_grid_chl[:,bad_depth]*np.nan
    mesh_grid_sig[:,bad_depth] = mesh_grid_sig[:,bad_depth]*np.nan
    mesh_grid_turb[:,bad_depth] = mesh_grid_turb[:,bad_depth]*np.nan
    mesh_grid_stats[:,bad_depth] = mesh_grid_stats[:,bad_depth]*np.nan

extent = (date_time.min(), date_time.max(), press_grid.max(), press_grid.min()) # extent of the plots


if args.image:

    fig = plt.figure()
    """
    ax = plt.subplot2grid((9,100), (0, 0), colspan=98)
    plt.plot(date2num(tw_time,'days since 1-1-1'), tw_wind,'k')
    ax.annotate('Wind Speed (m/s)', xy=(0, 1), xycoords='axes fraction', fontsize=6,
                xytext=(1, -5), textcoords='offset points',
                ha='left', va='top')
    #ax.set_ylabel('Wind Speed m/s', color='k',rotation=0, labelpad=-90,horizontalalignment='right')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_xlim([extent[0],extent[1]])
    ax.set_yticks(range(0,15,4))
    ax2 = ax.twinx()
    plt.plot(date2num(tw_time,'days since 1-1-1'), tw_td,'g')
    plt.plot(date2num(tw_time,'days since 1-1-1'), tw_temp,'r')
    #ax2.set_ylabel('Air Temp. (C)', color='r',rotation=0, labelpad=90,horizontalalignment='right')
    ax2.annotate('Air Temp. (C)', xy=(1, 1), xycoords='axes fraction', fontsize=6, color='r',
                xytext=(-1, -5), textcoords='offset points',
                ha='right', va='top')
    ax2.annotate('Dewpoint (C)', xy=(0.5, 1), xycoords='axes fraction', fontsize=6, color='g',
                xytext=(-0, -5), textcoords='offset points',
                ha='right', va='top')
    ax2.set_yticks(range(0,18,4))
    ax2.set_xlim([extent[0],extent[1]])
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    """
    ax = plt.subplot2grid((9,100), (1, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_t), extent=extent, cmap=cmocean.cm.thermal, vmin=-2.0, vmax=10.0, aspect='auto', alpha=0.85)
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Temperature (C)',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (2, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_s), extent=extent, cmap=cmocean.cm.haline, vmin=31.30, vmax=32.6, aspect='auto')
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Salinitity (PSU)',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (3, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_sig), extent=extent, cmap=cmocean.cm.dense, vmin=23.6, vmax=26.0, aspect='auto', alpha=0.85)
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('SigmaT (kg/m^3)',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (4, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_o), extent=extent, cmap=cmocean.cm.oxy, vmin=150, vmax=350, aspect='auto')
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Dissolved O2 Conc',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (5, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_osat), extent=extent, cmap=cmocean.cm.delta_r, vmin=60, vmax=120, aspect='auto')
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Dissolved O2 (% Sat)',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (6, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_chl), extent=extent, cmap=cmocean.cm.algae, vmin=0, vmax=15, aspect='auto')
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Chlor. A (ug/ml)',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (7, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_turb), extent=extent, cmap=cmocean.cm.turbid, vmin=0, vmax=50, aspect='auto')
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Turbidity (NTU)',rotation=0, labelpad=90,horizontalalignment='right')
    ax = plt.subplot2grid((9,100), (8, 0), colspan=100)
    cs = plt.imshow(np.transpose(mesh_grid_stats), extent=extent, cmap=cmocean.cm.gray_r, vmin=0, vmax=5, aspect='auto')
    cs.cmap.set_under('w')
    ax.xaxis.set_major_locator(DayLocator(bymonthday=15))
    ax.xaxis.set_minor_locator(DayLocator(bymonthday=[5,10,15,20,25,30]))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Points Per {0}m bin'.format(interval),rotation=0, labelpad=90,horizontalalignment='right')

    ax.xaxis.set_minor_formatter(DateFormatter('%d'))
    ax.xaxis.set_major_formatter(DateFormatter('%b %y'))
    ax.xaxis.set_tick_params(which='major', pad=15)

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*0.5, DefaultSize[1]*1.25) )
    plt.savefig('2018_BSITAE_prawler_image.svg', bbox_inches='tight', dpi = (300))
    plt.savefig('2018_BSITAE_prawler_image.png', bbox_inches='tight', dpi = (300))
    plt.close()

if args.csv_out in ['temperature','T_20','temp']:
    #at specified depth
    '''
    print "time,depth,T_20"
    for i in range(0,np.shape(mesh_grid_t)[0],1):
        datetemp = num2date(date_time[i]+1.,'Days since 0001-1-1')
        print "{0},{1},{2}".format(datetemp.strftime('%Y-%m-%d %H:%M:%S'),press_grid[70],mesh_grid_t[i,70])
    '''

    #all values as depthxtime array
    np.savetxt('SP03_temp.csv',mesh_grid_t,fmt='%.3f')
