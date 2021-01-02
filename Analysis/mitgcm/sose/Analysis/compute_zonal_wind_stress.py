import numpy as np
import xarray as xr



root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
print(project_name)
outname = project_name + '_zonal_mean_stress.nc'
fnames = root_folder+project_name+'/*oceTAUX*'
grname = root_folder+project_name+'/grid.nc'
list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)
gr = xr.open_dataset(grname)

taux = df.mean(('time','i_g'))
taux.to_netcdf(outname)
