import numpy as np 

gr = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')

# root_folder = '/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/rebuilt/'
root_folder = '/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_bdy7_geom3_rlx/rebuilt/'
root_folder = '/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_nemo4_2_06/rebuilt/'

ls1 = sorted(glob.glob(root_folder+'*grid_T.nc'))

ds = xr.open_dataset('BlackSea_mask.nc')
df = xr.open_dataset(ls1[0])
df = df.fillna(0)
mask = xr.where(df.thetao[0,:,:,:]!=0,1,0)
mask = mask*ds.mask
area = gr.e1t*gr.e2t
df = xr.open_mfdataset(ls1)

# temperature
# volvar = df.thetao*gr.e3t_0.isel(t=0).rename({'z': 'deptht'})*area.isel(t=0)
# vol = mask*gr.e3t_0.isel(t=0).rename({'z': 'deptht'})*area.isel(t=0)
# meanvar = volvar.sum(('x','y'))/vol.sum(('x','y'))
# temp1D = meanvar.resample(time_counter='1MS').mean('time_counter')
# temp1D.load()
# salinity

# df2 = df.so.where(df.so!=0)
volvar = df.so*mask*gr.e3t_0.isel(t=0).rename({'z': 'deptht'})*area.isel(t=0)
vol = mask*gr.e3t_0.isel(t=0).rename({'z': 'deptht'})*area.isel(t=0)
meanvar = volvar.sum(('x','y'))/vol.sum(('x','y'))
meanvar2 = volvar.sum(('x','y','deptht'))/vol.sum(('x','y','deptht'))
salt1D = meanvar.resample(time_counter='1MS').mean('time_counter')
salt1D.load()

bins = np.array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160,
                 180, 200, 300, 400, 500, 1000, 2000])
# fig = plt.figure(dpi=120)
# ax = plt.axes()
# plt.pcolormesh(np.arange(1,62), temp1D.deptht, np.transpose(np.copy(temp1D-temp1D[0,:])),
#                cmap='needJet2', vmin=-2, vmax=2)
# plt.semilogy()
# ax.invert_yaxis()
# ax.set_yticks(bins)
# ax.set_yticklabels(bins);

fig = plt.figure(dpi=120)
ax = plt.axes()
plt.pcolormesh(np.arange(1,salt1D.time_counter.shape[0]+1), salt1D.deptht, np.transpose(np.copy(salt1D-salt1D[0,:])),
               cmap='needJet2', vmin=-.5, vmax=.5)
plt.semilogy()
ax.invert_yaxis()
ax.set_yticks(bins)
ax.set_yticklabels(bins);


