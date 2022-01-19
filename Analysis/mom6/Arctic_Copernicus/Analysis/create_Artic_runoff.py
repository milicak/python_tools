import os
from matplotlib.path import Path as mpPath
import glob
from mpl_toolkits.basemap import Basemap

root_folder = '~/dataset/MOM6/Arctic_Copernicus/'

df = xr.open_dataset('~/dataset/CORE2/NIAF/runoff.daitren.iaf.20120419.nc')
gr = xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')
m = Basemap(projection='npstere',boundinglat=10,lon_0=0,resolution='l');
m.fillcontinents(color='coral',lake_color='aqua');
lon = np.copy(df.xc)
lat = np.copy(df.yc)
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

fout = root_folder + 'runoff.daitren.iaf.20120419_Arctic.nc'
df['runoff'] = df['runoff']*inside
df.to_netcdf(fout)
# add units kg/s/m^2 to the netcdf file
cmnds = 'ncatted -O -a units,runoff,c,c,kg/s/m^2 ' + fout
os.sys(cmnds)
# add Fill value uo the netcdf file
cmnds = 'ncatted -O -a _FillValue,runoff,c,d,-9999. ' + fout
os.sys(cmnds)


# ls1 = sorted(glob.glob('/cluster/shared/noresm/inputdata/ocn/jra55/IAF/JRA.v1.1.runoff_liquid*'))
# for filename in ls1:
#     print(filename[70:74])
#     ds = xr.open_dataset(filename)
#     ds['runoff'] = ds['runoff']*inside
#     fout = root_folder + 'JRA_Arctic_Copernicus_runoff_' + filename[70:74] + '.nc'
#     ds.to_netcdf(fout)
#
#
