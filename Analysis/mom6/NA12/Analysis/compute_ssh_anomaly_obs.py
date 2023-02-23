import numpy as np

ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/SSHUGVG/global/*.nc'))

df = xr.open_mfdataset(ls1)

sla_std = df.sla.std('time')
adt_std = df.adt.std('time')

ds = sla_std.to_dataset(name='std')
ds['atd'] = adt_std
ds.to_netcdf('/Volumes/A1/workdir/milicak/datasets/SSHUGVG/global_ssh_anomaly.nc')
