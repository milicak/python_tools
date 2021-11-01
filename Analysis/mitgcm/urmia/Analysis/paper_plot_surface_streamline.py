import xmitgcm
import cmocean 



ds = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/',
                             prefix='UVEL',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/',
                             prefix='VVEL',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

ds2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/',
                              prefix='UVEL',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/',
                              prefix='VVEL',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

ds3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/',
                              prefix='UVEL',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/',
                              prefix='VVEL',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

# Uvel1 = np.copy(ds.UVEL[99,0,:,:])
# Uvel2 = np.copy(ds2.UVEL[99,0,:,:])
# Vvel1 = np.copy(dt.VVEL[99,0,:,:])
# Vvel2 = np.copy(dt2.VVEL[99,0,:,:])

ds = ds.sel(time=slice('2019-04-01','2019-04-30'))
ds2 = ds2.sel(time=slice('2019-04-01','2019-04-30'))
ds3 = ds3.sel(time=slice('2019-04-01','2019-04-30'))
dt = dt.sel(time=slice('2019-04-01','2019-04-30'))
dt2 = dt2.sel(time=slice('2019-04-01','2019-04-30'))
dt3 = dt3.sel(time=slice('2019-04-01','2019-04-30'))

Uvel1 = np.copy(ds.UVEL[:,0,:,:].mean('time'))
Uvel2 = np.copy(ds2.UVEL[:,0,:,:].mean('time'))
Uvel3 = np.copy(ds3.UVEL[:,0,:,:].mean('time'))
Vvel1 = np.copy(dt.VVEL[:,0,:,:].mean('time'))
Vvel2 = np.copy(dt2.VVEL[:,0,:,:].mean('time'))
Vvel3 = np.copy(dt3.VVEL[:,0,:,:].mean('time'))

sp1 = np.sqrt(Uvel1**2+Vvel1**2) 
sp2 = np.sqrt(Uvel2**2+Vvel2**2) 
sp3 = np.sqrt(Uvel3**2+Vvel3**2) 

current_cmap = plt.cm.get_cmap(cmocean.cm.speed)
# current_cmap = plt.cm.get_cmap('viridis')
current_cmap.set_bad(color='gray')    
lon,lat = np.meshgrid(np.copy(ds.XC),np.copy(ds.YC)) 

fig, axes = plt.subplots(figsize=(21,6))
ax1 = plt.subplot2grid(shape=(1,3),loc=(0,0), colspan=1)
ax2 = plt.subplot2grid(shape=(1,3),loc=(0,1), colspan=1)
ax3 = plt.subplot2grid(shape=(1,3),loc=(0,2), colspan=1)
plt.tight_layout()


im1 = ax1.pcolormesh(ds.XC,ds.YC,ma.masked_equal(sp1,0),cmap=current_cmap,vmin=0,vmax=0.5)
fig.colorbar(im1, ax=ax1)
ax1.streamplot(lon, lat, Uvel1, Vvel1, color='k', density=2)   
im2 = ax2.pcolormesh(ds.XC,ds.YC,ma.masked_equal(sp3,0),cmap=current_cmap,vmin=0,vmax=0.5)
fig.colorbar(im2, ax=ax2)
ax2.streamplot(lon, lat, Uvel3, Vvel3, color='k', density=2)   
im3 = ax3.pcolormesh(ds.XC,ds.YC,ma.masked_equal(sp2,0),cmap=current_cmap,vmin=0,vmax=0.5)
fig.colorbar(im3, ax=ax3)
ax3.streamplot(lon, lat, Uvel2, Vvel2, color='k', density=2)   
ax1.set_xlabel('Longitude', fontsize=14);
ax1.set_ylabel('Latitude', fontsize=14);
ax2.set_xlabel('Longitude', fontsize=14);
ax2.set_ylabel('Latitude', fontsize=14);
ax3.set_xlabel('Longitude', fontsize=14);
ax3.set_ylabel('Latitude', fontsize=14);

ax1.text(45,37.15,'a)',fontsize=14) 
ax2.text(45,37.15,'b)',fontsize=14) 
ax3.text(45,37.15,'c)',fontsize=14) 

plt.savefig('paperfigs/urmia_surface_speed.png', bbox_inches='tight',format='png',dpi=300)


