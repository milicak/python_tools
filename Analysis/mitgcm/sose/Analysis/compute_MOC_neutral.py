import os
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
project_name = 'Exp01_0'
fname = root_folder + project_name + '/grid.nc'
gr  = xr.open_dataset(fname)

fname = root_folder + project_name + '/' + 'THETA_' + '2011_12_30.nc'
dft = xr.open_dataset(fname)
fname = root_folder + project_name + '/' + 'SALT_' + '2011_12_30.nc'
dfs = xr.open_dataset(fname)
fname = root_folder + project_name + '/' + 'VVELMASS_' + '2011_12_30.nc'
dfv = xr.open_dataset(fname)

dft = dft.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
dfs = dfs.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
dfv = dfv.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})

dfs1 = dfs.isel(time=0)
dft1 = dft.isel(time=0)
dfv1 = dfv.isel(time=0)
sigma2 = xr.apply_ufunc(gsw.sigma2, dfs1.SALT, dft1.THETA,
                            dask='parallelized', output_dtypes=[dfs1.SALT.dtype])
df1 = sigma2.to_dataset(name='sigma2')
df1['VVEL'] = df1.sigma2
df1['drF'] = df1.Z

lon = df1['XC'].values
lat = df1['YC'].values
zlev  = df1['Z'].values
me = xr.DataArray(np.copy(dfv1.VVELMASS), coords={'YC': lat, 'XC': lon,
                                'Z': zlev},
             dims=['Z', 'YC', 'XC'])
df1['VVEL'] = me
me = xr.DataArray(np.copy(gr.drF), coords={'Z': zlev},
             dims=['Z'])
df1['drF'] = me

# select bins for sigma2 or neutral
bins = np.arange(34, 38, 0.025)

df_rebinned = vertical_rebin_wrapper(df1,
                                     'sigma2',
                                     bins,
                                     dz_name='drF',
                                     vert_dim='Z')
df_rebinned = df_rebinned.fillna(0)
me = xr.DataArray(np.copy(gr.dxC), coords={'YC': lat, 'XC': lon},
             dims=['YC', 'XC'])
df_rebinned['dxC'] = me
voltr = df_rebinned.VVEL*df_rebinned.drF*df_rebinned.dxC





