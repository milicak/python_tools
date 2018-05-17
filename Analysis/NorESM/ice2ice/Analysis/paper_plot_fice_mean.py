import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from netCDF4 import num2date
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
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


def running_mean(x, N):
     cumsum = np.cumsum(np.insert(x, 0, 0))
     return (cumsum[N:] - cumsum[:-N]) / float(N)


#fyear = 1701
#lyear = 1922
fyear = 3301
lyear = 3600

gridfile = '/tos-project1/NS2345K/noresm/inputdata/ocn/micom/tnx1v1/20120120/grid.nc'
lon = xr.open_dataset(gridfile, decode_times=False)['plon']
lat = xr.open_dataset(gridfile, decode_times=False)['plat']
lon = np.copy(lon.data)
lat = np.copy(lat.data)
lon = np.copy(lon[:-1,:])
lat = np.copy(lat[:-1,:])

# 1 for atlantic_arctic_ocean region
# 2 for indian_pacific_ocean region
# 3 for global_ocean
chunks = (96,144)
xr_chunks = {'lat': chunks[-2], 'lon': chunks[-1]}
fname = '/tos-project1/NS4659K/milicak/data/Pacific2_fice.nc'
ctlname = '/tos-project1/NS4659K/milicak/data/ctl_fice.nc'
#data = xr.open_mfdataset(fname,chunks=xr_chunks)['TREFHT']
#dataref = xr.open_mfdataset(ctlname,chunks=xr_chunks)['TREFHT']
data = xr.open_dataset(fname, decode_times=False)['fice']
dataref = xr.open_dataset(ctlname, decode_times=False)['fice']
taux = np.copy(data.mean('time'))
taux = np.copy(taux[:-1,:])
tauxc = np.copy(dataref.mean('time'))
tauxc = np.copy(tauxc[:-1,:])
#
[lon1,lat1,ficereal] = enable_global(lon,lat,tauxc)
[lon,lat,ficediff] = enable_global(lon,lat,taux-tauxc)
ficediff = ficediff/100
#
plt.figure(figsize=(8,4))
m=Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=90,projection='cyl')
#m = Basemap(width=12000000,height=8000000,
#            resolution='l',projection='stere',\
#            lat_ts=40,lat_0=90,lon_0=0.)
#m = Basemap(projection='stere',boundinglat=60,lon_0=0,resolution='l')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(np.transpose(lon-110),np.transpose(lat)
                   ,np.transpose(np.ma.masked_invalid(ficereal))/100
                   ,shading='flat',cmap='jet',vmin=0,vmax=1,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="15%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
plt.tight_layout()
plt.savefig('Pacific2_fice_real.png',dpi=300,bbox_inches='tight')

plt.clf()
plt.figure(figsize=(8,4))
m=Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=90,projection='cyl')
#m = Basemap(width=12000000,height=8000000,
#            resolution='l',projection='stere',\
#            lat_ts=40,lat_0=90,lon_0=0.)
#m = Basemap(projection='stere',boundinglat=60,lon_0=0,resolution='l')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(np.transpose(lon-110),np.transpose(lat)
                   ,np.transpose(np.ma.masked_invalid(ficediff))
                   ,shading='flat',cmap='RdYlBu_r',vmin=-.5,vmax=.5,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="15%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
plt.tight_layout()
plt.savefig('Pacific2_fice_mean2.png',dpi=300,bbox_inches='tight')

plt.clf();
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
#
##gl.xlabels_top = False
##     ...: gl.ylabels_left = False
##     ...: gl.xlines = False
##     ...: gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
##     ...: gl.xformatter = LONGITUDE_FORMATTER
##     ...: gl.yformatter = LATITUDE_FORMATTER
##     ...: gl.xlabel_style = {'size': 15, 'color': 'gray'}
##     ...: gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
#
plt.pcolormesh(lon-110,lat,ficediff,vmin=-.25,vmax=.25,cmap='RdYlBu_r',
               shading='flat',transform=ccrs.PlateCarree())
plt.colorbar(ax=ax, shrink=.75)
plt.tight_layout()
plt.savefig('Pacific2_fice_mean.png',dpi=300,bbox_inches='tight')
#plt.close()
#
