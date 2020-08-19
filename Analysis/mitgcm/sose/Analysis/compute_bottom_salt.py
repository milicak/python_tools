import xmitgcm



ds = df.SALT.mean('time')
#indice
ind = xr.where(df.hFacC==0,0,1)
ind = np.int16(np.sum(ind,0))

sbottom = np.zeros((ind.shape[0],ind.shape[1]))

# x = xr.DataArray(your_i_indices, dims='points')
# y = xr.DataArray(your_j_indices, dims='points')
# z = xr.DataArray(ind, dims='points')
# your_slice = your_variable.isel(x=x,y=y)

for iy in np.arange(0,ind.shape[0]):
    for ix in np.arange(0,ind.shape[1]):
        sbottom[iy,ix]=ds.isel(Z=ind[iy,ix]-1,XC=ix,YC=iy)
        # sbottom[iy,ix] = np.copy(ds[ind[iy,ix]-1,iy,ix])


root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
project_name = 'Exp01_0s'
fnames = root_folder+project_name+'/'

dtime = 100 # 100 or 80

df = xmitgcm.open_mdsdataset(fnames,
                             prefix='SALT',geometry='sphericalpolar',
                            ref_date='2006-12-31 12:0:0',delta_t=dtime)

df = df.sel(time=slice('2013-01-01','2013-12-31'))

ds2 = df.SALT.mean('time')

ds4 = xr.open_dataset('Salt_2013.nc')
sbottom2 = np.zeros((ind.shape[0],ind.shape[1]))
for iy in np.arange(0,ind.shape[0]):
    for ix in np.arange(0,ind.shape[1]):
        sbottom2[iy,ix] = ds4.SALT.isel(Z=ind[iy,ix]-1,XC=ix,YC=iy)



lon, lat = np.meshgrid(np.copy(df.XC),np.copy(df.YC))
m = Basemap(projection='spstere',boundinglat=-60,lon_0=180,resolution='l');
m.fillcontinents(color='grey');
m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sbottom),cmap='jet',latlon=True);plt.colorbar();
plt.clim(34.4,34.8)



list = sorted(glob.glob('/archive2/milicak/mitgcm/sose/Exp01_0/*SALT*'))
df = xr.open_mfdataset(list)
df = df.sel(time=slice('2011-01-01','2011-12-31'))
ds = df.SALT.mean('time')
#indice
gr = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/grid.nc')
ind = xr.where(gr.hFacC==0,0,1)
ind = np.int16(np.sum(ind,0))
sbottom3 = np.zeros((ind.shape[0],ind.shape[1]))
for iy in np.arange(0,ind.shape[0]):
    for ix in np.arange(0,ind.shape[1]):
        sbottom2[iy,ix]=ds.isel(k=ind[iy,ix]-1,i=ix,j=iy)





