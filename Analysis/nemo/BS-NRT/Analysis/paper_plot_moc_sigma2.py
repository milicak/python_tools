import numpy as np

df = xr.open_dataset('moc_meridional_sigma2_BSe3r1_mean.nc')    
moc = np.copy(df.vol_sigma_tr)  
moc2 = np.reshape(moc,(215,50,2))
moc2 = moc2.sum(axis=2)
sigma2_bin = 0.5*(np.copy(df.sigma2_bin[1::2])+np.copy(df.sigma2_bin[::2])) 


moc2[moc2==0] = np.nan
cmap = plt.cm.get_cmap('RdBu_r') 
cmap.set_bad(color='grey')
plt.pcolormesh(df.lat,sigma2_bin,np.transpose(ma.masked_equal(moc2,0))*1e-6,cmap=cmap);plt.colorbar(); 
plt.clim(-0.1,0.1)  
