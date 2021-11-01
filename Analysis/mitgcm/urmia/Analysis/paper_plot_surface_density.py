import xmitgcm

def find_rho(ds,dt):
    rhoNil = 1237
    sBeta = 9.5e-4
    tAlpha = 3.5e-4
    dRho = 1237 - 1000
    refSalt = 278
    refTemp = 23
    rho = rhoNil*(sBeta*(ds.SALT-refSalt)-tAlpha*(dt.THETA-refTemp) ) + dRho
    return rho

ds = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/', prefix='SALT',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/', prefix='THETA',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

ds2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/', prefix='SALT',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/', prefix='THETA',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

ds3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/', prefix='SALT',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/', prefix='THETA',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

rho = find_rho(ds,dt)
rho2 = find_rho(ds2,dt2)
rho3 = find_rho(ds3,dt3)
rho = rho.where(rho>0,0)
rho2 = rho2.where(rho2>0,0)
rho3 = rho3.where(rho3>0,0)

current_cmap = plt.cm.get_cmap('YlOrRd')     
current_cmap.set_bad(color='gray')    

fig, axes = plt.subplots(figsize=(21,6))
ax1 = plt.subplot2grid(shape=(1,3),loc=(0,0), colspan=1)
ax2 = plt.subplot2grid(shape=(1,3),loc=(0,1), colspan=1)
ax3 = plt.subplot2grid(shape=(1,3),loc=(0,2), colspan=1)
plt.tight_layout()


im1 = ax1.pcolormesh(rho.XC,rho.YC,ma.masked_equal(rho.data[99,0,:,:]+20,rho[0,0,0,0]+20),cmap=current_cmap,vmin=26,vmax=200)
fig.colorbar(im1, ax=ax1)
im2 = ax2.pcolormesh(rho.XC,rho.YC,ma.masked_equal(rho3.data[99,0,:,:]+20,rho3[0,0,0,0]+20),cmap=current_cmap,vmin=26,vmax=200)
fig.colorbar(im2, ax=ax2)
im3 = ax3.pcolormesh(rho.XC,rho.YC,ma.masked_equal(rho2.data[99,0,:,:]+20,rho2[0,0,0,0]+20),cmap=current_cmap,vmin=26,vmax=200)
fig.colorbar(im3, ax=ax3)
ax1.set_xlabel('Longitude', fontsize=14);
ax1.set_ylabel('Latitude', fontsize=14);
ax2.set_xlabel('Longitude', fontsize=14);
ax2.set_ylabel('Latitude', fontsize=14);
ax3.set_xlabel('Longitude', fontsize=14);
ax3.set_ylabel('Latitude', fontsize=14);

ax1.text(45,37.15,'a)',fontsize=14) 
ax2.text(45,37.15,'b)',fontsize=14) 
ax3.text(45,37.15,'c)',fontsize=14) 

plt.savefig('paperfigs/urmia_surface_density.png', bbox_inches='tight',format='png',dpi=300)

