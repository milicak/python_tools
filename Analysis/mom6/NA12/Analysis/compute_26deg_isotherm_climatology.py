from xgcm import Grid
from mpl_toolkits.basemap import Basemap


dft = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/WOA13/woa13_decav_t00_04v2.nc',decode_times=False)
Temp_cr = 26.0
isotherm_depth = np.zeros((dft.lat.shape[0],dft.lon.shape[0]))
tmp = dft.t_an.where(dft.t_an>Temp_cr)
bb = xr.where(dft.t_an>Temp_cr,1,0)
itr = bb.sum('depth')
isotherm_depth = dft.depth[itr]
# create dataset                                             
dfs = xr.Dataset({
    'isotherm_depth': xr.DataArray(
                data   = isotherm_depth,
                dims   = ['time','lat','lon'],
                coords = {'time': dft.time,'lat': dft.lat, 'lon': dft.lon},
                attrs  = {
                   'units'     : 'meter'
                   }
               ),
           },
   )

dfs.to_netcdf('Isotherm_depth_climatology.nc')
