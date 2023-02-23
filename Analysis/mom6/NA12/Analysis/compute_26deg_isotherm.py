from xgcm import Grid
from mpl_toolkits.basemap import Basemap

gr = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_geometry.nc')
ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/*ocean_month_z*'))

# df.temp[0,:,500,800].interp(z_l=ds.zl,method="linear",kwargs={"fill_value": "extrapolate"})

Temp_cr = 26.0
# skip 1996 and 1997
ls1 = ls1[24:]
isotherm_depth = np.zeros((len(ls1),gr.geolon.shape[0],gr.geolon.shape[1]))
for ind,fname in enumerate(ls1):
    print(ind)
    df = xr.open_dataset(fname)
    tmp = df.thetao.where(df.thetao>Temp_cr)
    bb = xr.where(df.thetao>Temp_cr,1,0)
    itr = bb.sum('z_l')
    isotherm_depth[ind,:,:] = df.z_l[itr]

# dataset
time = pd.date_range("1998-01-01", freq="M", periods=len(ls1))
# create dataset                                             
dfs = xr.Dataset({
    'isotherm_depth': xr.DataArray(
                data   = isotherm_depth,
                dims   = ['time','lat','lon'],
                coords = {'time': time},
                attrs  = {
                   'units'     : 'meter'
                   }
               ),
           },
   )

dfs.to_netcdf('Isotherm_depth_MOM6.nc')
