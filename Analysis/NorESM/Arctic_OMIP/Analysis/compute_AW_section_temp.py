import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import os.path
from scipy.interpolate import interp1d
import ESMF
from netCDF4 import Dataset
import geopy.distance

# if __name__ == '__main__':
#     from dask.distributed import Client, LocalCluster
#     cluster = LocalCluster(n_workers=4, memory_limit=25e9
#                            ,local_dir='/cluster/work/users/milicak/OMIPs_diag/')
#     client = Client(cluster)


def create_src_grid(gridfile, data):
    '''

    '''
    grid_center_lon = np.copy(data.lon).flatten()
    grid_center_lat = np.copy(data.lat).flatten()
    grid_corner_lon = np.copy(data.lon_bnds)
    grid_corner_lon = grid_corner_lon.reshape(grid_center_lon.shape[0],4)
    grid_corner_lat = np.copy(data.lat_bnds)
    grid_corner_lat = grid_corner_lat.reshape(grid_center_lon.shape[0],4)
    mask=np.copy(data.thetao[0,0,:,:].fillna(0))
    mask[mask!=0]=1
    dataset = Dataset(gridfile, 'w', format='NETCDF4_CLASSIC')

    grid_size = dataset.createDimension('grid_size', grid_center_lon.shape[0])
    grid_corners = dataset.createDimension('grid_corners', 4)
    grid_rank = dataset.createDimension('grid_rank', 2)

    tmp = dataset.createVariable('grid_dims', 'i',('grid_rank',))
    tmp.long_name = 'grid_dims'
    tmp[:] = (df.lon.shape[1],df.lon.shape[0])

    tmp = dataset.createVariable('grid_imask', 'i',('grid_size',))
    tmp.units = 'unitless'
    tmp.long_name = 'grid_imask'
    tmp[:] = mask

    tmp = dataset.createVariable('grid_center_lon', 'd',('grid_size',))
    tmp.units = 'degrees'
    tmp.long_name = 'grid_center_lon'
    tmp[:] = grid_center_lon

    tmp = dataset.createVariable('grid_center_lat', 'd',('grid_size',))
    tmp.units = 'degrees'
    tmp.long_name = 'grid_center_lat'
    tmp[:] = grid_center_lat

    tmp = dataset.createVariable('grid_corner_lon', 'd',('grid_corners','grid_size'))
    tmp.units = 'degrees'
    tmp.long_name = 'grid_corner_lon'
    tmp[:] = grid_corner_lon

    tmp = dataset.createVariable('grid_corner_lat', 'd',('grid_corners','grid_size'))
    tmp.units = 'degrees'
    tmp.long_name = 'grid_corner_lat'
    tmp[:] = grid_corner_lat

    dataset.close()
    return


# section locations
lon_s4=np.array([17.6, 16.5, 16.05, 15.6, 15.1, 14.1, 13.0, 12.0, 10.0, 8.0, 4.0
             , 4.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0
             , 110.0, 120.0, 130.0, 140.0]);
lat_s4=np.array([69.0, 70.6, 71.3, 72.02, 72.8, 73.8, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0
    , 81.0, 81.8, 81.8, 82.6, 83.0, 83.2, 83.1, 82.8, 82.5, 81.8, 79.7, 78.2
    , 78.7, 79.7]);

x = np.linspace(1, lon_s4.shape[0], num=lon_s4.shape[0], endpoint=True)
f = interp1d(x,lon_s4)
g = interp1d(x,lat_s4)
xnew = np.linspace(1, lon_s4.shape[0], num=5*lon_s4.shape[0],
                   endpoint=True)

lon_s4new = f(xnew)
lat_s4new = g(xnew)

dnm = np.vstack([lat_s4new,lon_s4new])
dist = np.zeros(lon_s4new.shape)
for ind in range(0,lon_s4new.shape[0]-1):
    dist[ind+1] = dist[ind]+geopy.distance.distance(dnm[:,ind+1],dnm[:,ind]).km


coord_sys=ESMF.CoordSys.SPH_DEG
domask=True
# create locstream
locstream = ESMF.LocStream(lon_s4new.shape[0], name="Atlantic Inflow Section", coord_sys=coord_sys)
# appoint the section locations
locstream["ESMF:Lon"] = lon_s4new
locstream["ESMF:Lat"] = lat_s4new
if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_s4new.shape[0]), dtype=np.int32)


root_folder = '/tos-project1/NS9252K/CMIP6/'
omipn = 'omip1'
model = 'GFDL-CM4'
scn = 'r1i1p1f1'

# fnames = root_folder+omipn+'/'+ model+'/'+scn+'/thetao_Omon_'+model+'_'+omipn+'_'+scn+'_gn*'
# fnames = root_folder+omipn+'/'+ model+'/'+scn+'/so_Omon_'+model+'_'+omipn+'_'+scn+'_gn*'

fnames = root_folder+omipn+'/'+model+'/'+scn+'/thetao_Omon_'+model+'_'+omipn+'_'+scn+'_gn_'+'198801-200712.nc'
# fnames = root_folder+omipn+'/'+model+'/'+scn+'/so_Omon_'+model+'_'+omipn+'_'+scn+'_gn_'+'198801-200712.nc'

list=sorted(glob.glob(fnames))

CHUNKS = {'lev': 1}
# CHUNKS = {'y': 270, 'x': 360}



df = xr.open_mfdataset(list,chunks=CHUNKS)

dz = np.copy(df.lev)

# df = xr.open_mfdataset(list,chunks={'time': 20})
# df = xr.open_mfdataset(list)
# df = xr.open_mfdataset(list, parallel=True)

# create grid file
gridfile = model+'_ESMF_grd.nc'
if not os.path.isfile(gridfile):
    create_src_grid(gridfile, df)


# time mean of the thetao variable
# tmp = df.groupby('time.year').mean('time')
# kind = 0
# tmp=df.thetao[0,:,:,:].load()


secfield = np.zeros((lon_s4new.shape[0],df.lev.shape[0]))

for kind in range(0,df.lev.shape[0]):
    print('indice = ', kind)
    tmp = df.thetao[:,kind,:,:]
    tmp = tmp.groupby('time.year').mean('time')
    tmp = tmp.mean('year')
    # Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
    srcgrid = ESMF.Grid(filename=gridfile,
                     filetype=ESMF.FileFormat.SCRIP)

    srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
    # temp = np.copy(tmp[kind,:,:])
    # srcfield.data[:] = np.transpose(temp)
    srcfield.data[:] = np.transpose(tmp.data.compute())

    # create a field on the locstream
    dstfield = ESMF.Field(locstream, name='dstfield')
    dstfield.data[:] = 0.0

    # create an object to regrid data from the source to the destination field
    dst_mask_values=None
    if domask:
            dst_mask_values=np.array([0])

    regrid = ESMF.Regrid(srcfield, dstfield,
                        #regrid_method=ESMF.RegridMethod.NEAREST_STOD,
                        regrid_method=ESMF.RegridMethod.BILINEAR,
                        unmapped_action=ESMF.UnmappedAction.IGNORE,
                        dst_mask_values=dst_mask_values)

    # do the regridding from source to destination field
    dstfield = regrid(srcfield, dstfield)
    secfield[:,kind] = dstfield.data



secfield
del tmp
tmp = xr.DataArray(secfield, coords=[dist, dz], dims=['distance', 'dz'])
AW_section = tmp.to_dataset(name='temp_sec')
svname = model+'_'+omipn+'_'+scn+'_AW_section_temp.nc'
AW_section.to_netcdf(svname)
