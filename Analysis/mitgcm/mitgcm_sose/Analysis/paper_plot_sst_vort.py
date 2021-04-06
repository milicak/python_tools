import numpy as np

df = xr.open_dataset('SST_vort_SSS_Jan_2nd.nc')
lonc = np.copy(df.XC)
latc = np.copy(df.YC)
lonc, latc = np.meshgrid(lonc,latc)

#
plt.figure(figsize=(10,10))
m = Basemap(projection='spstere',boundinglat=-25,lon_0=0,resolution='l')
# lon, lat = m(lonc,latc);
m.bluemarble();
m.pcolormesh(lonc,latc,
             ma.masked_where(df.Depth==0,df.THETA[-1,:,:]),
             cmap='nice_gfdl',vmin=-2,vmax=28, latlon=True);
savename = 'gifs/sst_Jan2nd_snapshot.png'
plt.savefig(savename,bbox_inches='tight',format='png',dpi=150)



plt.figure(figsize=(10,10))
m = Basemap(projection='spstere',boundinglat=-25,lon_0=0,resolution='l')
m.bluemarble();
m.pcolormesh(lonc,latc,
             ma.masked_where(df.Depth==0,df.momVort3[-1,:,:]),
             cmap='RdBu_r',vmin=-5e-5,vmax=5e-5, latlon=True);
savename = 'gifs/vort_Jan2nd_snapshot.png'
plt.savefig(savename,bbox_inches='tight',format='png',dpi=150)
