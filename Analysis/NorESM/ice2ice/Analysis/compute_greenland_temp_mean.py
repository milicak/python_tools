import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import num2date
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
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
fyear = 3451
lyear = 3500
plt.figure()
root_folder = '/cluster/work/users/milicak/archive/'
#root_folder = '/tos-project1/NS4659K/chuncheng/cases_fram/'
root_folderref = '/tos-project1/NS4659K/chuncheng/cases_ice2ice/'
#expid = 'NBF1850_f19_tn11_test_mis3b_fwf3b_fram'
#expid = 'NBF1850_f19_tn11_test_mis3b_fwf3b_MI'
expidref = 'NBF1850_f19_tn11_test_mis3b_mixing3'
expid = 'NBF1850_f19_tn11_test_mis3b_mixing3_SO'
foldername = root_folder + expid + '/atm/hist/'
foldernameref = root_folderref + expidref + '/atm/hist/'
sdate="%c%4.4d%c" % ('*',fyear,'*')
freq = '*h0*'
list=sorted(glob.glob(foldername+freq+sdate))
listref=sorted(glob.glob(foldernameref+freq+sdate))
for year in xrange(fyear+1,lyear+1):
    sdate="%c%4.4d%c" % ('*',year,'*')
    list.extend(sorted(glob.glob(foldername+freq+sdate)))
    listref.extend(sorted(glob.glob(foldernameref+freq+sdate)))


chunks = (96,144)
xr_chunks = {'lat': chunks[-2], 'lon': chunks[-1]}
data = xr.open_mfdataset(list,chunks=xr_chunks)['TS']
dataref = xr.open_mfdataset(listref,chunks=xr_chunks)['TS']

dnm = data.mean('time')
dnmref = dataref.mean('time')

plt.clf();
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(data.lon,data.lat,dnm-dnmref,vmin=-5,vmax=5,cmap='RdYlBu_r',transform=ccrs.PlateCarree())
plt.colorbar(ax=ax, shrink=.62)
ax.coastlines()

