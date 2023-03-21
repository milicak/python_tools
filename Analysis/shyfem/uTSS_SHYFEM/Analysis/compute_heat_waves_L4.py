import numpy as np
from datetime import date
# Load marineHeatWaves definition module
import marineHeatWaves as mhw

df = xr.open_dataset('/work/opa/mi19918/Projects/nemo/Marmara/Marmara_global_L4_SST_2004_2020.nc')
df2 = xr.open_dataset('/work/opa/mi19918/Projects/nemo/Marmara/Marmara_global_L4_SST_2021.nc')
ds = xr.merge((df,df2))

# time for L4 satellite dataset
t = np.arange(date(2004,1,1).toordinal(),date(2021,12,31).toordinal()+1)
dates = [date.fromordinal(tt.astype(int)) for tt in t]

lons=np.array([27.525,28.825,28.825,27.525,27.525])
lats=np.array([40.525,40.525,40.925,40.925,40.525])


xstr = 20
xend = 46
ystr = 10
yend = 18

# 18 years times total points
count = np.zeros(((xend-xstr)*(yend-ystr),18))
duration = np.zeros(((xend-xstr)*(yend-ystr),18))
total_days = np.zeros(((xend-xstr)*(yend-ystr),18))
ind = 0
for jind in range(10,18):
    for iind in range(20,46):
        print(ind)
        sst = np.copy(ds.analysed_sst[:,jind,iind])
        # mhws, clim = mhw.detect(t, sst, climatologyPeriod=[2004,2010])
        mhws, clim = mhw.detect(t, sst, climatologyPeriod=[2017,2021])
        # mhws, clim = mhw.detect(t, sst)
        mhwBlock = mhw.blockAverage(t, mhws)
        count[ind] = mhwBlock['count']
        duration[ind] = mhwBlock['duration']
        total_days[ind] = mhwBlock['total_days']
        ind += 1



duration[np.isnan(duration)]=0

dnm=duration.mean(axis=1)
lon=ds.lon[xstr:xend]
lat=ds.lat[ystr:yend]
durationmean=dnm.reshape(8,26)
import cartopy.feature
import cartopy.crs as ccrs

ax = plt.axes(projection=ccrs.Mercator(central_longitude=180))
axpos = ax.get_position()
pos_x = axpos.x0+axpos.width + 0.01# + 0.25*axpos.width
pos_y = axpos.y0
cax_width = 0.04
cax_height = axpos.height
# Draw coastlines so we know where we are:
ax.coastlines()
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([26.5, 30, 40, 41.4], crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND,color='gray')
im=ax.pcolormesh(lon,lat,durationmean,transform=ccrs.PlateCarree())
plt.colorbar(im,ax=ax, shrink=.60)
