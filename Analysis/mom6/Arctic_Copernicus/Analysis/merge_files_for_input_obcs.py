'''This is code is written using
   https://github.com/ESMG/PyCNAL_regridding/blob/master/examples/SODA3.3.1/Creating_Initial_and_Boundary_conditions_from_SODA3.py'''
import os                                                  
import numpy as np                                         
import xesmf as xe                                         
import xarray as xr                                        
import scipy.io                                            
from scipy.io import savemat                               
from scipy.io import loadmat                               
# from mpl_toolkits.basemap import Basemap, shiftgrid        
import matplotlib.colors as colors                         
from scipy.signal import medfilt2d                         
import netCDF4                                             
import matplotlib.pyplot as plt                            
from scipy.interpolate import griddata                     
from matplotlib.path import Path                           
#for interpolation                                         
from scipy.spatial import cKDTree                          
from HCtFlood.kara import flood_kara                       
from PyCNAL_regridding import * 
                                                           
def velocity_at_corners(ds_u,ds_v):
    x=ds_u.x[-1].data+1;y=ds_v.y[-1].data+1
    #upper-right q points
    u_q=0.5*(ds_u.data+ds_u.data.roll(roll_coords='yh',yh=-1)).isel(xq=slice(1,x))
    #upper-right q points
    v_q=0.5*(ds_v.v+ds_v.v.roll(roll_coords='xh',xh=-1)).isel(yq=slice(1,y))
    ds_uvq = xr.Dataset({'u':u_q,'v':v_q},coords={'time':ds_u.time_counter,'lon':parent_grid['q'].x,'lat':parent_grid['q'].y,'angle_dx':parent_grid['q'].angle_dx})
    return ds_uvq


# dfu=xr.open_dataset('/data/products/OMIP/OMIP2_ORCA025/OMIP2.025d.01_1m_19580101_19581231_grid_U.nc')
# dfv=xr.open_dataset('/data/products/OMIP/OMIP2_ORCA025/OMIP2.025d.01_1m_19580101_19581231_grid_V.nc')
# ds_u=dfu['uo'];ds_v=dfv['vo']
# model_data['ds_uv']=velocity_at_corners(ds_u,ds_v)


# compression options
comp = dict(zlib=True, complevel=5)
root_folder = '/data/products/OMIP/OMIP2_ORCA025/'
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'

ls1 = sorted(glob.glob(root_folder+'*grid_T*'))
xstr = 305
xend = 366

for ind in range(xstr,xend):
    fnamet = ls1[ind]
    fnameu = fnamet[0:74]+'U.nc'
    fnamev = fnamet[0:74]+'V.nc'
    dft = xr.open_dataset(fnamet)
    dfu = xr.open_dataset(fnameu)
    dfv = xr.open_dataset(fnamev)
    dfu = dfu.rename({'nav_lon':'nav_lon_grid_T','nav_lat':'nav_lat_grid_T','x':'x_grid_T','y':'y_grid_T','depthu':'deptht'})
    dfv = dfv.rename({'nav_lon':'nav_lon_grid_T','nav_lat':'nav_lat_grid_T','x':'x_grid_T','y':'y_grid_T','depthv':'deptht'})
    ds = dft[['thetao','so','zos']]
    ds['uo'] = dfu.uo
    ds['vo'] = dfv.vo
    outfile = mom_dir + fnamet[34:]
    encoding = {var: comp for var in ds.data_vars}
    print(outfile)
    ds.to_netcdf(outfile, encoding=encoding)

# compute angle 
# Approximate angles using centered differences in interior                     
# lon = np.copy(df.nav_lon_grid_T)
# lat = np.copy(df.nav_lat_grid_T)
# angle = np.zeros((df.nav_lon_grid_T.shape[0],df.nav_lon_grid_T.shape[1]))
# angle[:,1:-1] = np.arctan( (lat[:,2:]-lat[:,:-2]) /                             
#                           ((lon[:,2:]-lon[:,:-2])*np.cos(np.deg2rad(lat[:,1:-1]))) )
# # Approximate angles using side differences on left/right edges                 
# angle[:,0] = np.arctan( (lat[:,1]-lat[:,0]) / ((lon[:,1]-lon[:,0])*np.cos(np.deg2rad(lat[:,0]))) )
# angle[:,-1] = np.arctan( (lat[:,-1]-lat[:,-2]) /                                
#                         ((lon[:,-1]-lon[:,-2])*np.cos(np.deg2rad(lat[:,-1]))) ) 
#
# ds_u=dfu['uo'];ds_v=dfv['vo']
# model_data['ds_uv']=velocity_at_corners(ds_u,ds_v)
