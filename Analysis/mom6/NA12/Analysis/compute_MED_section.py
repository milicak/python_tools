import numpy as np
import glob
from mpl_toolkits.basemap import Basemap


root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'

gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')

ls1 = sorted(glob.glob(root_folder + '*ocean_month_z*'))
# ls1 = sorted(glob.glob(root_folder + '*annual.nc*'))

# Take the last 10 years
# ls1 = ls1[-10:]
ls1 = ls1[8:]
df = xr.open_mfdataset(ls1)

ds1 = df.thetao.isel(yh=887).mean('time')
ds2 = df.so.isel(yh=887).mean('time')

ds = ds1.to_dataset(name='temp')
ds['salt'] = ds2
ds.to_netcdf('MED_section_38_5_MOM6.nc')

