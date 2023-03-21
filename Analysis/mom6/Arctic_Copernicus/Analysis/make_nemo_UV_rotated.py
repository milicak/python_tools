import numpy as np
import os                                                  
import numpy as np                                         
import xesmf as xe                                         
import glob
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

root_folder = '/data/products/OMIP/OMIP2_ORCA025/'
fname1 = 'OMIP2.025d.01_1m_19580101_19581231_grid_T.nc'
df = xr.open_dataset(root_folder + fname1)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'                     
ls1 = sorted(glob.glob(root_folder+'*grid_T*'))
xstr = 305
xend = 366

# compute angle 
# Approximate angles using centered differences in interior                     
lon = np.copy(df.nav_lon_grid_T)
lon1 = np.copy(df.nav_lon_grid_T)
lat = np.copy(df.nav_lat_grid_T)
cos_lat = np.cos(np.radians(lat))
angle1 = np.zeros(lat.shape)
angle2 = np.zeros(lat.shape)
angle1[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
angle1[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
angle1[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
lon = np.where(lon < 0., lon+360, lon)
angle2[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
angle2[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
angle2[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
angle = np.maximum(angle1, angle2)
lon = lon1

ds2 = xr.Dataset({
            "angle": (["y_grid_T", "x_grid_T"], angle),
        },
        coords={
            "nav_lon_grid_T": (["y_grid_T", "x_grid_T"],
                               np.copy(df.nav_lon_grid_T)),
            "nav_lat_grid_T": (["y_grid_T", "x_grid_T"],
                               np.copy(df.nav_lat_grid_T)),
        },
    )


# compression options
comp = dict(zlib=True, complevel=5)
for ind0 in range(xstr,xend):
    fname1 = ls1[ind0][-44:]
    df = xr.open_dataset(mom_dir + fname1)
    # apply rotation
    ue=np.cos(angle)*df.uo-np.sin(angle)*df.vo
    vn=np.sin(angle)*df.uo+np.cos(angle)*df.vo
    df['ue'] = ue
    df['vn'] = vn
    outname = mom_dir + fname1[:-3]+'_rotated.nc'
    df = df.drop_vars(('thetao','so','uo','vo','zos'))
    encoding = {var: comp for var in df.data_vars}
    print(outname)
    df.to_netcdf(outname, encoding=encoding)






