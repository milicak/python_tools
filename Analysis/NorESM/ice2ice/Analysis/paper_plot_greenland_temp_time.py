import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import num2date
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import glob

plt.ion()

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
              '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
              'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
              'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
              'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30,
                                      31, 30, 31],
              'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
              '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
              '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
              (year % 100 == 0) and (year % 400 != 0) and
              (year < 1583)):
            leap = False
    return leap


def get_dpm(time, calendar='standard'):
        """
            return a array of days per month corresponding to the months
        provided in `months`
            """
        month_length = np.zeros(len(time), dtype=np.int)
        cal_days = dpm[calendar]
        print cal_days
        print time.month

        for i, (month, year) in enumerate(zip(time.month, time.year)):
            month_length[i] = cal_days[month]
            if leap_year(year, calendar=calendar):
                month_length[i] += 1
        return month_length


def running_mean(x, N):
     cumsum = np.cumsum(np.insert(x, 0, 0))
     return (cumsum[N:] - cumsum[:-N]) / float(N)


#fyear = 1701
#lyear = 1922
fyear = 3301
lyear = 3600


# 1 for atlantic_arctic_ocean region
# 2 for indian_pacific_ocean region
# 3 for global_ocean
chunks = (96,144)
xr_chunks = {'lat': chunks[-2], 'lon': chunks[-1]}
fname = '/tos-project1/NS4659K/milicak/data/Pacific2_airtemp.nc'
ctlname = '/tos-project1/NS4659K/milicak/data/ctl_airtemp.nc'
#data = xr.open_mfdataset(fname,chunks=xr_chunks)['TREFHT']
#dataref = xr.open_mfdataset(ctlname,chunks=xr_chunks)['TREFHT']
data = xr.open_dataset(fname, decode_times=False)['TREFHT']
dataref = xr.open_dataset(ctlname, decode_times=False)['TREFHT']
lat1 = 65.0;
lat2 = 80.0;
lon1 = -55.0+360;
lon2 = -15.0+360;
jnd1 = np.int(np.min(np.where(data.lat>=lat1)));
jnd2 = np.int(np.max(np.where(data.lat<=lat2)));
ind1 = np.int(np.min(np.where(data.lon>=lon1)));
ind2 = np.int(np.max(np.where(data.lon<=lon2)));
#
area = xr.open_dataset('/tos-project1/NS4659K/milicak/data/area_cam.nc',
                       decode_times=False)['area']
areatime = np.tile(np.copy(area.data),((lyear-fyear+1)*12,1,1))
#
sat = {}
satannual = {}
dnm = np.nansum(data[:,jnd1:jnd2+1,ind1:ind2+1]*areatime[:,jnd1:jnd2+1,ind1:ind2+1],axis=(1,2))
dnmarea = np.nansum(areatime[:,jnd1:jnd2+1,ind1:ind2+1],axis=(1,2))
sat['SO'] = dnm/dnmarea
dnm = np.nansum(dataref[:,jnd1:jnd2+1,ind1:ind2+1]*areatime[:,jnd1:jnd2+1,ind1:ind2+1],axis=(1,2))
sat['ctrl'] = dnm/dnmarea
dnm = np.copy(sat['SO'])
dnm = np.reshape(dnm,(lyear-fyear+1,12))
satannual['SO'] = np.mean(dnm,axis=1)
dnm = np.copy(sat['ctrl'])
dnm = np.reshape(dnm,(lyear-fyear+1,12))
satannual['ctrl'] = np.mean(dnm,axis=1)
satrunmean = running_mean(satannual['SO']-satannual['ctrl'],10)
#
#cmap = plt.get_cmap('Set1')
#plt.figure(figsize=(8,4))
#line1, = plt.plot(satannual['SO']-satannual['ctrl'],'r',label='Air Temp diff')
#line2, = plt.plot(satrunmean,'k',label='10 year mean')
#plt.ylabel('Temperature [$^\circ$C]')
#plt.xlabel('Time [years]')
#plt.legend(loc='lower right')
#plt.tight_layout()
#
satrunmean = running_mean(satannual['SO']-satannual['ctrl'],5)
plt.figure(figsize=(8,4))
line1, = plt.plot(satannual['SO']-satannual['ctrl'],'r',label='Air Temp diff')
line2, = plt.plot(satrunmean,'k',label='5 year mean')
plt.legend(loc='lower right')
plt.ylabel('Temperature [$^\circ$C]')
plt.xlabel('Time [years]')
#plt.tight_layout()
plt.savefig('Pacific2_temp_greenland.png',dpi=300,bbox_inches='tight')

