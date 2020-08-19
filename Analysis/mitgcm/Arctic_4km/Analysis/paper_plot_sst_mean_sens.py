import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
import matplotlib as mpllib

# AtlorPac = 'Atl'
AtlorPac = 'Atlv2'
# AtlorPac = 'Pac'

df = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_0_SST.nc')
if AtlorPac == 'Atl':
    dfw = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_1_SST.nc')
elif AtlorPac == 'Atlv2':
    dfw = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_3_SST.nc')
else:
    dfw = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_2_SST.nc')



lon = np.copy(df.XC)
lat = np.copy(df.YC)

fig = plt.figure(figsize=(8,24))
# fig = plt.figure()
# fig.canvas.draw()
plt.tight_layout()
m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l');
x,y = m(-45, 34.5)

# 1-5 years
ax1 = plt.subplot2grid((6,4),(0,1),colspan = 2, rowspan = 2)
dnm = df.sel(year=slice(1992,1996))
sstc = dnm.mean('year')
dnm = dfw.sel(year=slice(1992,1996))
sstw = dnm.mean('year')
m.ax = ax1
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sstw.sst-sstc.sst),
                   cmap='RdBu_r',vmin=-2,
           vmax=2,latlon=True, rasterized=True)
cb = m.colorbar(im1,"right", size="5%", pad="4%")
plt.text(x,y,'a)', fontsize=14);

# 6-10 years
ax2 = plt.subplot2grid((6,4),(2,0),colspan = 2, rowspan = 2)
dnm = df.sel(year=slice(1997,2001))
sstc = dnm.mean('year')
dnm = dfw.sel(year=slice(1997,2001))
sstw = dnm.mean('year')
m.ax = ax2
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sstw.sst-sstc.sst),
                   cmap='RdBu_r',vmin=-2,
           vmax=2,latlon=True, rasterized=True)
plt.text(x,y,'b)', fontsize=14);
# 11-15 years
ax3 = plt.subplot2grid((6,4),(2,2),colspan = 2, rowspan = 2)
dnm = df.sel(year=slice(2002,2006))
sstc = dnm.mean('year')
dnm = dfw.sel(year=slice(2002,2006))
sstw = dnm.mean('year')
m.ax = ax3
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sstw.sst-sstc.sst),
                   cmap='RdBu_r',vmin=-2,
           vmax=2,latlon=True, rasterized=True)
plt.text(x,y,'c)', fontsize=14);
# 16-20 years
ax4 = plt.subplot2grid((6,4),(4,0),colspan = 2, rowspan = 2)
dnm = df.sel(year=slice(2007,2011))
sstc = dnm.mean('year')
dnm = dfw.sel(year=slice(2007,2011))
sstw = dnm.mean('year')
m.ax = ax4
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sstw.sst-sstc.sst),
                   cmap='RdBu_r',vmin=-2,
           vmax=2,latlon=True, rasterized=True)
plt.text(x,y,'d)', fontsize=14);
# 21-25 years
ax5 = plt.subplot2grid((6,4),(4,2),colspan = 2, rowspan = 2)
dnm = df.sel(year=slice(2012,2016))
sstc = dnm.mean('year')
dnm = dfw.sel(year=slice(2012,2016))
sstw = dnm.mean('year')
m.ax = ax5
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sstw.sst-sstc.sst),
                   cmap='RdBu_r',vmin=-2,
           vmax=2,latlon=True, rasterized=True)
plt.text(x,y,'e)', fontsize=14);
plt.tight_layout()
fig.subplots_adjust(wspace=-.5)


if AtlorPac == 'Atl':
    plt.savefig('paperfigs/sst_mean_Atlantic.png', bbox_inches='tight',format='png',dpi=300)
elif AtlorPac == 'Atlv2':
    plt.savefig('paperfigs/sst_mean_Atlanticv2.png', bbox_inches='tight',format='png',dpi=300)
else:
    plt.savefig('paperfigs/sst_mean_Pacific.png', bbox_inches='tight',format='png',dpi=300)


