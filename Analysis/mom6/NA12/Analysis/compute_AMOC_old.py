from xgcm import Grid
import xesmf as xe
from mpl_toolkits.basemap import Basemap
import os


root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'
deg_rad=np.pi/180.

gr = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_geometry.nc')
gr2 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/ocean_hgrid.nc')
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})
angle_dx = gr2.angle_dx[1::2,1::2]
out = xr.Dataset()
out['lon'] = gr.geolon
out['lat'] = gr.geolat
# out = gr.rename({'geolon': 'lon', 'geolat': 'lat'})

df0 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/19981001.ocean_month_z_1998_12.nc')
df0 = df0.merge(gr)

# 2D grid
grid = Grid(df0, coords={'X': {'center': 'xh', 'outer': 'xq'},
                        'Y': {'center': 'yh', 'outer': 'yq'},
                         }, periodic=[])

ls1 = sorted(glob.glob(root_folder + '*ocean_month_z*'))
# last 20 years
ls1 = ls1[24:]
ds = xr.Dataset()
ds['lon'] = gr.xh[1:-1:2]
ds['lat'] = np.arange(-27.5,80,0.5)
moc = np.zeros((len(ls1),df0.z_l.shape[0],ds.lat.shape[0]))

for idx,fname in enumerate(ls1):
    print(fname)
    df = xr.open_dataset(fname)
    df = df.merge(gr)
    # compute mass transport on rho points
    urho = grid.interp(df.umo, 'X', boundary='fill')
    vrho = grid.interp(df.vmo, 'Y', boundary='fill')
    # ue = np.cos(angle_dx.data*deg_rad)*urho.data-np.sin(angle_dx.data*deg_rad)*vrho.data
    vn = np.sin(angle_dx.data*deg_rad)*urho+np.cos(angle_dx.data*deg_rad)*vrho
    # verticalsum = np.flip(np.nancumsum(np.flip(vn,1),axis=1),1)
    # build regridder
    if os.path.exists('bilinear_1844x1678_215x838.nc'):
        regridder = xe.Regridder(out, ds, 'bilinear',reuse_weights=True,
                periodic=False)
    else:
        regridder = xe.Regridder(out, ds, 'bilinear',reuse_weights=False,
                periodic=False)
    dr_out = regridder(vn)
    verticalsum = dr_out.reindex(z_l=dr_out.z_l[::-1]).cumsum('z_l').reindex(z_l=dr_out.z_l[:])
    moc[idx,:,:] = -verticalsum.sum(('time','xh'))*1e-9
    # moc[idx,:,:] = -verticalsum.sum(axis=(0,3))*1e-9


meanmoc = moc.mean(axis=0)
time = pd.date_range("1998-01-01", freq="M", periods=len(ls1))
# create dataset                                             
dfs = xr.Dataset({                                           
    'AMOC': xr.DataArray(                             
                data   = moc,                              
                dims   = ['time','z_l','lat'],                          
                coords = {'z_l': np.copy(df.z_l),'time': time,'lat':
                    np.copy(ds.lat)},       
                attrs  = {                                   
                    'units'     : 'Sv'                      
                    }                                        
                ),                                           
            },                                               
    )                                                        

dfs.to_netcdf('AMOC_MOM6.nc')
