from xgcm import Grid
from mpl_toolkits.basemap import Basemap

gr = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_geometry.nc')
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})

# df0 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/20010101.ocean_daily.nc')
# df0 = df0.merge(gr)
ds = xr.open_dataset('MOM6_mean_velocities.nc')
ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/*ocean_daily*.nc'))
root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/uprime/'
# skip 1996 and 1997

for ind, fname in enumerate(ls1):
    print(ind)
    df = xr.open_dataset(fname)
    df = df.merge(gr)
    # 2D grid
    grid = Grid(df, coords={'X': {'center': 'xh', 'outer': 'xq'},
                            'Y': {'center': 'yh', 'outer': 'yq'},
                             }, periodic=[])

    urho = grid.interp(df.ssu, 'X', boundary='fill')
    vrho = grid.interp(df.ssv, 'Y', boundary='fill')
    upr = urho - ds.urho_mean 
    vpr = vrho - ds.vrho_mean 
    df1 = upr.to_dataset(name='uprime')
    df1['vprime'] = vpr
    fout = root_folder + ls1[ind][-23:-8] + 'uprime.nc'
    df1.to_netcdf(fout)
    df1.close()
    del df1


    # if ind==0:
    #     upr = urho - ds.urho_mean
    #     vpr = vrho - ds.vrho_mean
    # else:
    #     upr = xr.concat((upr,urho - ds.urho_mean),'time')
    #     vpr = xr.concat((vpr,vrho - ds.urho_mean),'time')
    #


# umean = urho.mean('time')
# vmean = vrho.mean('time')
# upr = urho - umean
# vpr = vrho - vmean
#
# EKE = upr**2+vpr**2
# EKE = EKE.mean('time')
#
# plt.pcolormesh(gr.geolon,gr.geolat,EKE);plt.colorbar()
# plt.axis([-85,-30,25,55])
# plt.clim(0,0.7)
#
# MKE = umean**2+vmean**2
# plt.pcolormesh(gr.geolon,gr.geolat,MKE);plt.colorbar()
# plt.axis([-85,-30,25,55])
# plt.clim(0,1.5)
#
# vorticity = ( - grid.diff(df0.ssu * df0.dxCu, 'Y', boundary='fill')
#               + grid.diff(df0.ssv * df0.dyCv, 'X', boundary='fill') ) / df0.Aq
