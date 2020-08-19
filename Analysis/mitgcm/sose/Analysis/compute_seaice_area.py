import numpy as np
import xmitgcm

root_folder = '/archive2/milicak/mitgcm/sose/'
expid = 'Exp01_0'
outname = expid + '_SIarea.nc'
grname = root_folder+expid+'/grid.nc'
gr = xr.open_dataset(grname)

fnames = root_folder + expid + '/*SIarea_*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)
dfn = df.SIarea.rename({'i': 'XC', 'j': 'YC'})
dnm = dfn.data*gr.rA
icearea = dnm.sum(axis=(1,2))
aa = icearea.compute()

dsn = xr.Dataset({'SIarea': ('time', aa), 'time': df.time})
dsn.to_netcdf(outname)


# time = pd.date_range('2007-01-01', freq='1D', periods=2563)
