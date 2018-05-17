import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
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


def enable_global(tlon,tlat,data):
  """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
  tlon=tlon+abs(ma.max(tlon)); tlon=tlon+360
  # stack grids side-by-side (in longitiudinal direction), so
  # any range of longitudes may be plotted on a world map.
  tlon = np.concatenate((tlon,tlon+360),1)
  tlat = np.concatenate((tlat,tlat),1)
  data = ma.concatenate((data,data),1)
  tlon = tlon-360.
  return tlon, tlat, data


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


gridfile = '/tos-project1/NS2345K/noresm/inputdata/ocn/micom/tnx1v1/20120120/grid.nc'
#fyear = 1701
#lyear = 1922
fyear = 3301
lyear = 3600
#root_folder = '/tos-project1/NS4659K/milicak/'
#root_folder = '/cluster/work/users/milicak/archive/'
#root_folder = '/tos-project1/NS4659K/chuncheng/cases_fram/'
root_folder = '/tos-project1/NS4659K/milicak/'
root_folderref = '/tos-project1/NS4659K/chuncheng/cases_ice2ice/'
#expid = 'NBF1850_f19_tn11_test_mis3b_fwf3b_fram'
#expid = 'NBF1850_f19_tn11_test_mis3b_fwf3b_MI'
expidref = 'NBF1850_f19_tn11_test_mis3b_mixing3'
expid = 'NBF1850_f19_tn11_test_mis3b_mixing3_Pacific2'
foldername = '/ocn/hist/'
foldername = root_folder + expid + '/ocn/hist/'
foldernameref = root_folderref + expidref + '/ocn/hist/'
sdate="%c%4.4d%c" % ('*',fyear,'*')
freq = '*hm*'
list=sorted(glob.glob(foldername+freq+sdate))
listref=sorted(glob.glob(foldernameref+freq+sdate))
for year in xrange(fyear+1,lyear+1):
    sdate="%c%4.4d%c" % ('*',year,'*')
    list.extend(sorted(glob.glob(foldername+freq+sdate)))
    listref.extend(sorted(glob.glob(foldernameref+freq+sdate)))


fig = plt.figure()
chunks = (190,180)
#chunks = (385,360)
xr_chunks = {'x': chunks[-1], 'y': chunks[-2]}
area = xr.open_mfdataset(gridfile,chunks=xr_chunks)['parea']
area = np.copy(area.data)
area = np.copy(area[:-1,:])
data = xr.open_mfdataset(list,chunks=xr_chunks)['fice']
dataref = xr.open_mfdataset(listref,chunks=xr_chunks)['fice']
fice = np.copy(data.data)
ficeref = np.copy(dataref.data)
#fice = np.copy(data.mean('time'))
#ficeref = np.copy(dataref.mean('time'))
fice = np.copy(fice[:,:-1,:])
ficeref = np.copy(ficeref[:,:-1,:])
ice_cr = 0.15;
fice = fice/100
ficeref = ficeref/100
# apply critical value
fice[np.where(fice<ice_cr)] = 0.0
fice[np.where(fice>=ice_cr)] = 1.0
ficeref[np.where(ficeref<ice_cr)] = 0.0
ficeref[np.where(ficeref>=ice_cr)] = 1.0

tt = fice.shape
area = np.tile(area,(tt[0],1, 1))

fice = fice*area
ficeref = ficeref*area

fice = fice[:,191:,:]
fice = np.nansum(fice,axis=(1,2))

ficeref = ficeref[:,191:,:]
ficeref = np.nansum(ficeref,axis=(1,2))


plt.plot(np.mean(fice.reshape((fice.size/12,12)),1)*1e-12,'b',label='SPG')
plt.plot(np.mean(ficeref.reshape((ficeref.size/12,12)),1)*1e-12,'k',label='CTRL')
plt.legend(loc='lower right')
