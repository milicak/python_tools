import numpy as np
from pandas.tseries.offsets import DateOffset



root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'
df = xr.open_dataset(root_folder + 'TS_woa_Arctic_Copernicus.nc',decode_times=False)
gridname = root_folder + 'ocean_geometry.nc'
gr = xr.open_dataset(gridname)
lon = np.copy(gr.geolon)
lat = np.copy(gr.geolat)

ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[12:24]
for itr in range(3,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])


ds = xr.open_mfdataset(ls2[24:],decode_times=False)
time = pd.date_range("1998-01-15", freq=DateOffset(months=1), periods=ds.Time.shape[0])
ds['Time'] = time

dfs = ds.thetao[:,10,:,:].mean('Time')

# ind = 10 for mom6 and 26 for WOA13

# m = Basemap(width=2600000,height=2300000,
#             resolution='l',projection='stere',
#             lat_ts=50,lat_0=70,lon_0=0.)
m = Basemap(projection='npstere',boundinglat=37,lon_0=0,resolution='i');
lon1,lat1 = m(-132,26)

fig, axes = plt.subplots(figsize=(9,6))
ax1 = plt.subplot2grid(shape=(1,4),loc=(0,0),colspan=2)
ax2 = plt.subplot2grid(shape=(1,4),loc=(0,2),colspan=2)
plt.tight_layout()

m.ax=ax1
m.drawcoastlines()
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(gr.D==0,dfs[:,:]),
                   cmap='nice_gfdl',vmin=-1,vmax=12,latlon=True, rasterized=True)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x0-0.05,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im1, cax=cbar_ax, ticklocation='left')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r'$^\circ$C',rotation=0,y=1.02,labelpad=-30,fontsize=14)
ax1.text(lon1,lat1,'a)',fontsize=14)

m.ax=ax2
m.drawcoastlines()
m.fillcontinents(color='grey');
im2 = m.pcolormesh(lon,lat,dfs-np.copy(df.temp_woa[26,:,:]),
                   cmap=plt.cm.get_cmap('RdBu_r',16),vmin=-3,vmax=3,latlon=True, rasterized=True)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$^\circ$C',rotation=0,y=1.07,labelpad=-40,fontsize=14)
ax2.text(lon1,lat1,'b)',fontsize=14);

plt.savefig('paperfigs/temp_250m_bias.png', bbox_inches='tight',format='png',dpi=300)


