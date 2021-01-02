import os
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr


root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
outname = project_name + '_MOC_depth.nc'
fname = root_folder + project_name + '/grid.nc'
gr  = xr.open_dataset(fname)

list = sorted(glob.glob(root_folder + project_name + '/' + 'VVELMASS_*'))

dfv = xr.open_mfdataset(list)

dfv = dfv.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})

voltrV = dfv.VVELMASS*gr.dxG*gr.drF


Trx = voltrV.sum(('XC'))
Trxmean = Trx.mean('time')
# reverse z coordinate for different cumsum
Trxmeanreverse_z = Trxmean.reindex(Z=Trxmean.Z[::-1])

Trxmeanreverse_zsum = Trxmeanreverse_z.cumsum('Z')
Trxmeanreverse_zsum = Trxmeanreverse_zsum.reindex(Z=Trxmeanreverse_zsum.Z[::-1])
Trxmeanreverse_zsum = -Trxmeanreverse_zsum

ds = Trxmeanreverse_zsum.to_dataset(name='MOC_depth')
ds.to_netcdf(outname)

