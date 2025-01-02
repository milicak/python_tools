import xarray as xr
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# gr = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/EastMed/OUT/ocean_geometry.nc')
# df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/EastMed/OUT/ocean_daily_v1.nc')

# to plot in altay
gr = xr.open_dataset('ocean_geometry.nc')
df = xr.open_dataset('ocean_daily_v1.nc')

ssh = np.copy(df.zos)
time = pd.date_range("2023-02-06-04-17", freq="30S", periods=359)

# Erdemli
lon1 = 34.32793 ; lat1 = 36.61112;
dist1 = abs(gr.geolon[0,:]-lon1)
xind = np.copy(dist1.argmin())
dist1 = abs(gr.geolat[:,0]-lat1)
yind = np.copy(dist1.argmin())
xind += 2
yind -= 2
ds = df.ssh[:,yind,xind]
ds = ds.to_dataset(name='Erdemli')

# Arsuz
lon1 = 35.88519 ; lat1 = 36.41559
dist1 = abs(gr.geolon[0,:]-lon1)
xind = np.copy(dist1.argmin())
dist1 = abs(gr.geolat[:,0]-lat1)
yind = np.copy(dist1.argmin())
xind -= 3
ds['Arsuz'] = df.ssh[:,yind,xind]

# Gazimagus
lon1 = 33.93092 ; lat1 = 35.17399
lon1 = 33.9501 ; lat1 = 35.1234
dist1 = abs(gr.geolon[0,:]-lon1)
xind = np.copy(dist1.argmin())
dist1 = abs(gr.geolat[:,0]-lat1)
yind = np.copy(dist1.argmin())
yind += 1
ds['Gazimagusa'] = df.ssh[:,yind,xind]

# Tasucu
lon1 = 33.83622 ; lat1 = 36.28146
dist1 = abs(gr.geolon[0,:]-lon1)
xind = np.copy(dist1.argmin())
dist1 = abs(gr.geolat[:,0]-lat1)
yind = np.copy(dist1.argmin())
yind -= 1
xind += 1
ds['Tasucu'] = df.ssh[:,yind,xind]
ds['time'] = time

ds.to_netcdf('ssh_stations.nc')

# Girne
# 33.32036 35.34170
# Bozyazi
# lon1 = 32.94131 ; lat1 = 36.09742


#Creating the map object
fig = plt.figure(figsize=(9,4))
plt.plot(time,ds.Erdemli)
plt.xlabel('Time [hours]')
plt.ylabel('ssh [meter]')
fname = 'paperfigs/Antakya_tsunami1_ssh_Erdemli.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)
