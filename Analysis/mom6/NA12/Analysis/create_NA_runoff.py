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

lon_pac = np.array([ 9048029.34206251, 12573138.58634935, 12065283.86471481,
       11348312.49299545, 11497681.52877031, 11467807.72161534,
       11228817.26437556, 11019700.61429074, 10750836.34989598,
       10571593.50696614, 10392350.6640363 ,  9944243.55671171,
        9675379.29231695,  9077903.14921748,  8390805.5846531 ,
        8420679.39180807, 9048029.34206251])

lat_pac = np.array([3283837.76860507, 3283837.76860507, 5375004.2694532 ,
       5882858.99108774, 6540082.74849715, 7167432.69875159,
       7376549.3488364 , 7316801.73452646, 7346675.54168143,
       7585665.99892121, 8003899.29909084, 8003899.29909084,
       8242889.75633062, 8392258.79210549, 6420587.51987726,
       5584120.91953801, 3283837.76860507])

vertices2 = np.transpose(np.array([lon_pac,lat_pac]))
path2 = mpPath(vertices2)
inside2 = path2.contains_points(pnts)
inside2 = np.reshape(inside2,((ny,nx)))
inside2 = np.double(inside2)
inside[inside2==1] = 0

# black sea
lon_bs = np.array([22967935.510735497, 24767590.534039382,
                   24755934.737256326,22972597.82944872,
                  22967935.510735497])
lat_bs = np.array([10458202.148715083, 10446546.351932026,
                    11211166.620900515,11157549.955698457,
                  10458202.148715083])

vertices3 = np.transpose(np.array([lon_bs,lat_bs]))
path3 = mpPath(vertices3)
inside3 = path3.contains_points(pnts)
inside3 = np.reshape(inside3,((ny,nx)))
inside3 = np.double(inside3)
inside[inside3==1] = 0

# East africa
lon_afr = np.array([23302593.001464494,24203086.558646895,
                    24192263.318776917,23389178.92042434,
                   23302593.001464494])

lat_afr = np.array([3379250.4559573424, 3372756.512035354,
                     4552489.657863263, 4580630.081525214,
                   3379250.4559573424])

vertices4 = np.transpose(np.array([lon_afr,lat_afr]))
path4 = mpPath(vertices4)
inside4 = path4.contains_points(pnts)
inside4 = np.reshape(inside4,((ny,nx)))
inside4 = np.double(inside4)
inside[inside4==1] = 0

fout = root_folder + 'runoff.daitren.iaf.20120419_NA12.nc'
fout = root_folder + 'runoff.daitren.clim.20120419_NA12.nc'
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