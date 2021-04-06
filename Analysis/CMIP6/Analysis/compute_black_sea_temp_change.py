import gcsfs



def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds

##################################################################################
#ssp585
###############

gra = xr.open_dataset('/okyanus/users/dcetin/area/areacello_Ofx_GFDL-CM4_piControl_r1i1p1f1_gr.nc')
area = gra.areacello



sst_raw3_hist = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Omon/thetao/gr/')
t3_hist = sst_raw3_hist.so[:,0,130:140,28:43]
t3_hist = t3_hist.groupby('time.year').mean('time')
t3_hist = t3_hist.sel(year=slice('1990','2010'))
mean_year3_hist = t3_hist.mean("year")
#mean_year3_hist=mean_year3_hist.to_dataset(name = 'sst')
#mean_year3_hist.to_netcdf('BS_1990_2010_hist_sst_meantime.nc')
BS_1990_2010_hist_sst_meantime = xr.open_dataset('BS_1990_2010_hist_sst_meantime.nc')
#plt.pcolormesh(BS_1990_2010_hist_sst_meantime.sst)
t3_hist = t3_hist* area[130:140,28:43]
t3_hist_area = area[130:140,28:43]
t3_hist = t3_hist.sum(['lon','lat'])/t3_hist_area.sum(['lon','lat'])
#t3_hist=t3_hist.to_dataset(name = 'sst')
#t3_hist.to_netcdf('BS_1990_2010_sst_hist.nc')
sst_BS_2040_2060_ssp585 = xr.open_dataset("BS_1990_2010_hist_sst.nc")


plt.figure(figsize=(12,8))
plt.plot(np.arange(2040,2061), BS_1990_2010_hist_sst.sst, color='red', linewidth=2, label='')
plt.title('Sea Surface Temperature of Black Sea Region for historical Simulation (NCAR/CESM2-WACCM)')
plt.xlabel('Year')
plt.ylabel('Temperature (C)')
plt.show()

fig=plt.figure(figsize=(13, 6))
m = Basemap(projection='cyl', resolution='h',
            llcrnrlat=39, urcrnrlat = 48,
            llcrnrlon=26, urcrnrlon = 45)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='lightgray')
#m.drawcountries(linewidth=1)
m.drawparallels(np.arange(-90., 120., 5.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180, 180., 5.), labels=[0, 0, 0, 1])
m.pcolormesh(BS_1990_2010_hist_sst_meantime.lon, BS_1990_2010_hist_sst_meantime.lat,
             BS_1990_2010_hist_sst_meantime.sst, shading='interp', latlon=True, cmap='RdBu_r')
cb = m.colorbar(size='3%', pad='2%')
cb.set_label('SST (\u2103)', fontsize=14)
mpl.rcParams.update({'font.size': 9})
plt.title('Sea Surface Temperature of Black Sea Region for 1990-2010 (NCAR/CESM2-WACCM/ssp585)', fontsize=14, fontweight='bold', pad=30);
from pylab import *
fontsize = 14
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    #tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
plt.xlabel('Longitude', fontsize=14, labelpad=30)
plt.ylabel('Latitude', fontsize=14, labelpad=30)
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2
#clim(15,20)



