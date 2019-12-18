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




root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'

project_name = 'Exp01_0'

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
fnames = root_folder+project_name+'/'+'*VVELMASS*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)

dfv = xr.open_dataset('/archive/milicak/dataset/SOSE/1_over_3_degree/bsose_i105_2008to2012_3day_Vvel.nc')
dfs = xr.open_dataset('/archive/milicak/dataset/SOSE/1_over_3_degree/bsose_i105_2008to2012_3day_Salt.nc')
dft = xr.open_dataset('/archive/milicak/dataset/SOSE/1_over_3_degree/bsose_i105_2008to2012_3day_Theta.nc')
dfg = xr.open_dataset('/archive/milicak/dataset/SOSE/1_over_3_degree/bsose_i105_2008to2012_3day_GAMMA.nc')
dfg.rename({'iTIME': 'time', 'iDEPTH': 'Z', 'iLAT': 'YG','iLON': 'XC',})
# df = dft.merge(dfs)
df = dfg.merge(dfv)

dfs1=dfs.isel(time=0)
dft1=dft.isel(time=0)
dfv1=dfv.isel(time=0)

bins = np.arange(24,28,0.2)
bins = np.arange(20, 40, 0.5)
df_rebinned = vertical_rebin_wrapper(dfv1,
                                     'sigma2',
                                     bins,
                                     dz_name='drF',
                                     vert_dim='Z')

sigma2 = xr.apply_ufunc(gsw.sigma2, dfs1.SALT, dft1.THETA,
                            dask='parallelized', output_dtypes=[dfs1.SALT.dtype])
df1=xr.DataArray(sigma2,name='sigma2')

# dfv1['sigma2']=sigma2
voltrV = dfv1.VVEL*dfv1.hFacS*dfv1.dxG*dfv1.drF

df_rebinned = vertical_rebin(voltrV, sigma2, bins, voltrV.drF,vert_dim='Z')


# plt.pcolormesh(df2.YG,df2.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()

