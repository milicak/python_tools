import numpy as np

ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/*ocean_daily*.nc'))

# skip 1996 and 1997
df = xr.open_mfdataset(ls1[8:])
sla_std = df.zos.std('time')
ds = sla_std.to_dataset(name='zosstd')

ds.to_netcdf('NA12_ssh_anomaly.nc')
