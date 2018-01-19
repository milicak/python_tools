import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import num2date
import matplotlib.pyplot as plt
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

fyear = 1701
lyear = 1950
root_folder = '/cluster/work/users/milicak/archive/'
expid = 'NBF1850_f19_tn11_test_mis3b_fwf3b_MI'
foldername = '/ocn/hist/'
foldername = root_folder + expid + '/ocn/hist/'
sdate="%c%4.4d%c" % ('*',fyear,'*')
freq = '*hy*'
list=(glob.glob(foldername+freq+sdate))
for year in xrange(fyear+1,lyear+1):
    sdate="%c%4.4d%c" % ('*',year,'*')
    list.extend(glob.glob(foldername+freq+sdate))


chunks = (360,385)
xr_chunks = {'x': chunks[-1], 'y': chunks[-2]}
data = xr.open_mfdataset(list,chunks=xr_chunks)['templvl']
data1 = xr.open_mfdataset(list,chunks=xr_chunks)['salnlvl']
zr = xr.open_mfdataset(list,chunks=xr_chunks)['depth']
temp = np.copy(data.data[:,:,315,107])
salt = np.copy(data1.data[:,:,315,107])

time = np.linspace(fyear,lyear,lyear-fyear+1)
plt.clf()
plt.pcolormesh(time,-zr,np.transpose(temp),cmap='gist_ncar');plt.colorbar()
plt.ylim((-1400,0))
plt.figure()
plt.pcolormesh(time,-zr,np.transpose(salt),cmap='gist_ncar');plt.colorbar()
plt.ylim((-1400,0))
