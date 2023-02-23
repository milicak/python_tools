from xgcm import Grid
from mpl_toolkits.basemap import Basemap

gr = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_geometry.nc')
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})

# df0 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/20010101.ocean_daily.nc')
# df0 = df0.merge(gr)
ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/*ocean_daily*.nc'))
# skip 1996 and 1997
df = xr.open_mfdataset(ls1[8:])

ds1 = df.ssu.mean('time')
ds1.load()
ds2 = df.ssv.mean('time')
ds2.load()
ds = ds1.to_dataset(name='ssu')
ds['ssv'] = ds2
ds = ds.merge(gr)

# 2D grid
grid = Grid(ds, coords={'X': {'center': 'xh', 'outer': 'xq'},
                        'Y': {'center': 'yh', 'outer': 'yq'},
                         }, periodic=[])

urho = grid.interp(ds.ssu, 'X', boundary='fill')
vrho = grid.interp(ds.ssv, 'Y', boundary='fill')

ds3 = urho.to_dataset(name='urho_mean')
ds3['vrho_mean'] = vrho

ds3.to_netcdf('MOM6_mean_velocities.nc')
