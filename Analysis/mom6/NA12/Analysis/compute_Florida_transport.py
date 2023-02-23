import numpy as np
import xarray as xr
import glob


root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'

ls1 = sorted(glob.glob(root_folder + '*Florida*' ))
df = xr.open_mfdataset(ls1)

# compute transport in Sv
Fl_tr = df.vmo.sum(('z_l', 'yq_sub01', 'xh_sub02'))*1e-9
Fl_trm = Fl_tr.resample(time='M').mean('time')

time = pd.date_range("1996-01-15", freq="M", periods=Fl_trm.shape[0])
ds = Fl_trm.to_dataset(name='Florida_tr')
ds['time'] = time
ds.to_netcdf('Florida_Strait_transport.nc')

# plt.plot(time,Fl_trm)
# plt.grid()
# plt.ylim(10,40)
# plt.xlabel('time')
# plt.ylabel('Florida Transport [Sv]')
