import os                                                  
import numpy as np                                         
import xesmf as xe                                         
import xarray as xr                                        
import scipy.io                                            
from scipy.io import savemat                               
from scipy.io import loadmat                               
from mpl_toolkits.basemap import Basemap, shiftgrid        
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

# create a mask file first
# mask = np.ones(df.thetao[0,:,:,:].shape) 
# mask[np.where(df.thetao[0,:,:,:] == 0.0)] = 0
# ds = xr.Dataset({
#             "mask": (["deptht", "y_grid_T", "x_grid_T"], mask),
#         },
#         coords={
#             "nav_lon_grid_T": (["y_grid_T", "x_grid_T"], df.nav_lon_grid_T),
#             "nav_lat_grid_T": (["y_grid_T", "x_grid_T"], df.nav_lat_grid_T),
#             "depth":(["deptht"], df.deptht ),
#         },
#     )
# ds.to_netcdf('nemo0_25degree_mask.nc')

# ---------- define a domain target on MOM grid ---------------------
# Nx and Ny should be the hgrid sizes so double the model grid
Nx=4000
Ny=3500
domain = obc_segment('domain',path_regional_grid,istart=0,iend=Nx,jstart=0,  jend=Ny)

# ---------- define variables on each segment ------------------
temp_domain = obc_variable(domain,'temp',geometry='surface',obctype='radiation')
salt_domain = obc_variable(domain,'salt',geometry='surface',obctype='radiation')
ssh_domain  = obc_variable(domain,'ssh' ,geometry='line'   ,obctype='flather')


maskfile = mom_dir + 'nemo0_25degree_mask.nc' 

interp_t2s = temp_domain.interpolate_from(root_folder + fname1,'thetao',frame=0,depthname='deptht',
                                               coord_names=['nav_lon_grid_T','nav_lat_grid_T'],method='bilinear',
                                               missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  

interp_s2s = salt_domain.interpolate_from(root_folder + fname1,'so',frame=0,depthname='deptht',
                             coord_names=['nav_lon_grid_T','nav_lat_grid_T'],missing_value=0,maskfile=maskfile,
                             maskvar='mask',autocrop=False,method='bilinear')

interp_ssh2s = ssh_domain.interpolate_from(root_folder + fname1,'zos',frame=0,
                             coord_names=['nav_lon_grid_T','nav_lat_grid_T'],missing_value=0,maskfile=maskfile,
                             maskvar='mask',autocrop=False,method='bilinear')

interp_ssh2s = ssh_domain.interpolate_from(root_folder + fname1,'zos',frame=0,
                                           coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                           autocrop=False,method='bilinear',missing_value=0) 


# ---------- list segments and variables to be written -------
list_segments = [domain]
list_variables = [ssh_domain,temp_domain,salt_domain]

#----------- time --------------------------------------------
time = timeobject(15.5)
time.units = 'days since 1958-01-01'
time.calendar = 'gregorian'

# ---------- write to file -----------------------------------
write_ic_file(list_segments,list_variables,[],time,output=mom_dir + 'NEMO_TS_IC.nc') 
# write_ic_file(list_segments,list_variables,list_vectvariablestime,output=mom_dir + 'NEMO_TS_IC.nc') 


