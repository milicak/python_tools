import numpy as np
import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
expid = 'Exp01_0'
outname = expid + '_SIarea.nc'

fnames = root_folder + expid + '/Southern_Ocean_ctrl_monthly_ocn_SIarea_*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)
dnm = df.SIarea*df.rA
icearea = dnm.sum(('XC','YC'))

ds = icearea.to_dataset(name='SIarea')
ds.to_netcdf(outname)


df = xmitgcm.open_mdsdataset('/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/Exp01_0s',prefix='SIarea',ref_date="2006-12-31 12:0:0",delta_t=100)
dnm = df.SIarea*df.rA
icearea = dnm.sum(('XC','YC'))
outname = 'Exp01_0s_SIarea.nc'
ds = icearea.to_dataset(name='SIarea')
ds.to_netcdf(outname)
