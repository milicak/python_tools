import os
import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'

from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import gsw
from xhistogram.xarray import histogram

def vertical_rebin(data, bin_data, bins, dz, vert_dim="st_ocean"):
    nanmask = np.isnan(data)
    # Should we also check the bin data for nans?
    full_sum = histogram(
        bin_data.where(~nanmask),
        bins=[bins],
        weights=(data * dz).where(~nanmask),
        dim=[vert_dim],
    )
    return full_sum

def vertical_rebin_wrapper(
    ds,
    bin_data_name,
    bins,
    dz_name="dz",
    vert_dim="st_ocean",
    return_average=True,
    debug=False,
):
    """A wrapper for the core functionality in `vertical_rebin`.
    Accepts datasets and calculates the average over the new depth coordinates.
    """
    ds = ds.copy()
    ds_rebinned = xr.Dataset()

    ones = xr.ones_like(ds[dz_name])

    dz_rebinned = vertical_rebin(
        ones,
        ds[bin_data_name],
        bins,
        ds[dz_name],
        vert_dim=vert_dim,
    )
    for var in ds.data_vars:
        ds_rebinned[var] = vertical_rebin(
            ds[var], ds[bin_data_name], bins, ds[dz_name], vert_dim=vert_dim
        )
    if return_average:
        ds_rebinned = (
            ds_rebinned / dz_rebinned
        )  # this might cause a lot of overhead...i can try to deactivate if the save fails.

    ds_rebinned[dz_name] = dz_rebinned

    return ds_rebinned



root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
fname = root_folder + project_name + '/grid.nc'
gr  = xr.open_dataset(fname)
# select bins for sigma2 or neutral
bins = np.arange(34, 38, 0.025)

list = sorted(glob.glob('/archive2/milicak/mitgcm/sose/Exp01_0/THETA*.nc'))
comp = dict(zlib=True, complevel=5)

for ind in np.arange(0,len(list)):
    inname = list[ind][-13:]
    fname = root_folder + project_name + '/' + 'SIGMA2_' + inname
    dfs = xr.open_dataset(fname)
    fname = root_folder + project_name + '/' + 'VVELMASS_' + inname
    dfv = xr.open_dataset(fname)
    dfv = dfv.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
    dfv1 = dfv.isel(time=0)
    dfs['VVEL'] = dfs.sigma2
    dfs['drF'] = dfs.Z
    lon = dfs['XC'].values
    lat = dfs['YC'].values
    zlev  = dfs['Z'].values
    tmp1 = xr.DataArray(np.copy(gr.drF), coords={'Z': zlev},
                 dims=['Z'])
    tmp2 = xr.DataArray(np.copy(gr.dxC), coords={'YC': lat, 'XC': lon},
             dims=['YC', 'XC'])
    me = xr.DataArray(np.copy(dfv1.VVELMASS), coords={'YC': lat, 'XC': lon,
                                    'Z': zlev},
                 dims=['Z', 'YC', 'XC'])
    dfs['VVEL'] = me
    dfs['drF'] = tmp1
    dfs['Z']=np.copy(gr.Z)
    df_rebinned = vertical_rebin_wrapper(dfs,
                                         'sigma2',
                                         bins,
                                         dz_name='drF',
                                         vert_dim='Z')
    df_rebinned = df_rebinned.fillna(0)
    df_rebinned['dxC'] = tmp2
    voltr = df_rebinned.VVEL*df_rebinned.drF*df_rebinned.dxC
    voltr = voltr.to_dataset(name='vol_sigma_tr')
    voltr = voltr.sum('XC')
    fname = root_folder + project_name + '/' + 'VVEL_SIGMA2_' + inname
    encoding = {var: comp for var in voltr.data_vars}
    print(fname)
    voltr.to_netcdf(fname, encoding=encoding)





