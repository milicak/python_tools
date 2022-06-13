import os
from matplotlib.path import Path as mpPath
import glob
from mpl_toolkits.basemap import Basemap

root_folder = '~/dataset/MOM6/NA12/'

df = xr.open_dataset('~/dataset/CORE2/NIAF/runoff.daitren.iaf.20120419.nc')


df = xr.open_dataset('~/dataset/CORE2/NYF_v2.0/runoff.daitren.clim.v2011.02.10.nc')

gr = xr.open_dataset('~/dataset/MOM6/NA12/ocean_hgrid.nc')

m = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.fillcontinents(color='coral',lake_color='aqua');

lon_baltic = np.array([21243565.22693396, 23519277.961993467,
                      23049686.44523516, 21869687.249278378,
                       21020810.27667682, 21243565.22693396])

lat_baltic = np.array([11317523.690144135, 11738951.974414412,
                      12310890.360209791, 12178441.470867703,
                      11937625.308427544, 11317523.690144135])

lon = np.copy(df.xc)
lat = np.copy(df.yc)
lon1, lat1 = m(lon,lat)
vertices = np.transpose(np.array([lon_baltic,lat_baltic]))
pnts = np.transpose(np.array([lon1.flatten(),lat1.flatten()]))
path = mpPath(vertices)
inside = path.contains_points(pnts)
ny = lon.shape[0]
nx = lon.shape[1]
inside = np.reshape(inside,((ny,nx)))
inside = np.double(inside)


fout = root_folder + 'runoff.daitren.iaf.20120419_Baltic.nc'
fout = root_folder + 'runoff.daitren.clim.20120419_Baltic.nc'
df['runoff'] = df['runoff']*inside
df.to_netcdf(fout)

# first remove fill value just in case
cmnds = 'ncatted -h -O -a _FillValue,,d,, ' + fout

# add units kg/s/m^2 to the netcdf file
cmnds = 'ncatted -O -a units,runoff,c,c,kg/s/m^2 ' + fout
os.sys(cmnds)
# add Fill value uo the netcdf file
cmnds = 'ncatted -O -a _FillValue,runoff,c,d,-9999. ' + fout
os.sys(cmnds)
cmnds = 'ncatted -O -a _FillValue,time,c,d,9.96920996838687e+36 ' + fout

# ls1 = sorted(glob.glob('/cluster/shared/noresm/inputdata/ocn/jra55/IAF/JRA.v1.1.runoff_liquid*'))
# for filename in ls1:
#     print(filename[70:74])
#     ds = xr.open_dataset(filename)
#     ds['runoff'] = ds['runoff']*inside
#     fout = root_folder + 'JRA_Arctic_Copernicus_runoff_' + filename[70:74] + '.nc'
#     ds.to_netcdf(fout)
#
#
