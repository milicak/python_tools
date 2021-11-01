import xmitgcm

ds = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/', prefix='S',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
ds2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/', prefix='S',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
ds3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/', prefix='S',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')



# current_cmap = plt.cm.get_cmap('needJet2')     
current_cmap = plt.cm.get_cmap('gist_earth')     
current_cmap.set_bad(color='gray')    
depth = np.copy(ma.masked_equal(ds.Depth[:,:],0))  
depth[depth==0] = np.nan   
depth2 = np.copy(ma.masked_equal(ds2.Depth[:,:],0))  
depth2[depth2==0] = np.nan   
depth3 = np.copy(ma.masked_equal(ds3.Depth[:,:],0))  
depth3[depth3==0] = np.nan   


fig, axes = plt.subplots(figsize=(21,6))
ax1 = plt.subplot2grid(shape=(1,3),loc=(0,0), colspan=1)
ax2 = plt.subplot2grid(shape=(1,3),loc=(0,1), colspan=1)
ax3 = plt.subplot2grid(shape=(1,3),loc=(0,2), colspan=1)
plt.tight_layout()

im1 = ax1.pcolormesh(ds.XC,ds.YC,depth,cmap=current_cmap,vmin=0,vmax=2)
fig.colorbar(im1, ax=ax1)
im2 = ax2.pcolormesh(ds.XC,ds.YC,depth3,cmap=current_cmap,vmin=0,vmax=2)
fig.colorbar(im2, ax=ax2)
im3 = ax3.pcolormesh(ds.XC,ds.YC,depth2,cmap=current_cmap,vmin=0,vmax=2)
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

plt.savefig('paperfigs/urmia_bathy.png', bbox_inches='tight',format='png',dpi=300)


