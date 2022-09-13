import os
from matplotlib.path import Path as mpPath
import glob
from mpl_toolkits.basemap import Basemap
import xarray

root_folder = '~/dataset/MOM6/Arctic_Copernicus/'
year1 = 1996
fname = '~/dataset/GLOFAS/glofas-era5_' + str(year1) + '.nc'
df = xr.open_dataset(fname)
gr = xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')


uparea = xarray.open_dataarray('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/upArea.nc')

# Find river end points by looking for local maxima in upstream area.
uparea = uparea.fillna(0).values
points = np.zeros_like(uparea)
window = 2  # look with +- this number of grid points
ni, nj = uparea.shape
for i in range(window, ni-window):
    for j in range(window, nj-window):
        sub = uparea[i-window:i+window+1, j-window:j+window+1]
        point = uparea[i, j]
        # A river end point has a reasonably large upstream area
        # and is a local maximum
        if point > 1e6 and sub.max() == point:
            points[i, j] = 1

m = Basemap(projection='npstere',boundinglat=10,lon_0=0,resolution='l');
# m.fillcontinents(color='coral',lake_color='aqua');

lon = np.copy(df.lon)
lat = np.copy(df.lat)
lon, lat = np.meshgrid(lon,lat)
lon1, lat1 = m(lon,lat)
lon_mom = np.concatenate((gr.x[-1,::-1],gr.x[::-1,0],gr.x[0,:],gr.x[:,-1]))
lat_mom = np.concatenate((gr.y[-1,::-1],gr.y[::-1,0],gr.y[0,:],gr.y[:,-1]))
lon_mom1, lat_mom1 = m(lon_mom,lat_mom)
vertices = np.transpose(np.array([lon_mom1.flatten(),lat_mom1.flatten()]))
pnts = np.transpose(np.array([lon1.flatten(),lat1.flatten()]))
path = mpPath(vertices)
inside = path.contains_points(pnts)
ny = lon.shape[0]
nx = lon.shape[1]
inside = np.reshape(inside,((ny,nx)))
inside = np.double(inside)


# Convert m3/s to kg/m2/s
# Borrowed from https://xgcm.readthedocs.io/en/latest/xgcm-examples/05_autogenerate.html
distance_1deg_equator = 111000.0
dlon = dlat = 0.1  # GloFAS grid spacing
dx = dlon * np.cos(np.deg2rad(df.lat)) * distance_1deg_equator
dy = ((df.lon * 0) + 1) * dlat * distance_1deg_equator
glofas_area = dx * dy


for year in range(1996,2018):
    print(year)
    files = [f'/okyanus/users/milicak/dataset/GLOFAS/glofas-era5_{year}.nc' for year
             in [year-1, year, year+1]]

    glofas = (
        xarray.open_mfdataset(files, combine='by_coords')
        .sel(time=slice(f'{year-1}-12-31 00:00:00', f'{year+1}-01-01 12:00:00'))
        .dis24.chunk({'time': -1})
    )
    glofas_kg = glofas * 1000.0 / glofas_area
    glofas_kg = glofas_kg*inside*points
    #glofas_kg = glofas_kg.rename({'dis24': 'runoff'})
    #glofas_kg['area'] = glofas_area
    #glofas_kg['runoff'] = glofas_kg['runoff']*inside*points
    fout = root_folder + 'glofas-era5_Arctic_Copernicus_' + str(year) + '.nc'
    glofas_kg=glofas_kg.to_dataset(name='runoff')
    glofas_kg['lon']=glofas_kg.lon-360
    glofas_kg.to_netcdf(fout)
    glofas_kg.close()
    # first remove fill value just in case
    cmnds = 'ncatted -h -O -a _FillValue,,d,, ' + fout
    os.system(cmnds)
    # add units kg/s/m^2 to the netcdf file
    cmnds = 'ncatted -O -a units,runoff,c,c,kg/s/m^2 ' + fout
    os.system(cmnds)
    # add Fill value uo the netcdf file
    cmnds = 'ncatted -O -a _FillValue,runoff,c,d,-9999. ' + fout
    os.system(cmnds)
    cmnds = 'ncatted -O -a _FillValue,time,c,d,9.96920996838687e+36 ' + fout
    os.system(cmnds)


# ls1 = sorted(glob.glob('/cluster/shared/noresm/inputdata/ocn/jra55/IAF/JRA.v1.1.runoff_liquid*'))
# for filename in ls1:
#     print(filename[70:74])
#     ds = xr.open_dataset(filename)
#     ds['runoff'] = ds['runoff']*inside
#     fout = root_folder + 'JRA_Arctic_Copernicus_runoff_' + filename[70:74] + '.nc'
#     ds.to_netcdf(fout)
#
#
