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


dfu=xr.open_dataset('/data/products/OMIP/OMIP2_ORCA025/OMIP2.025d.01_1m_19580101_19581231_grid_U.nc')
dfv=xr.open_dataset('/data/products/OMIP/OMIP2_ORCA025/OMIP2.025d.01_1m_19580101_19581231_grid_V.nc')
ds_u=dfu['uo'];ds_v=dfv['vo']
model_data['ds_uv']=velocity_at_corners(ds_u,ds_v)



root_folder = '/data/products/OMIP/OMIP2_ORCA025/'
fname1 = 'OMIP2.025d.01_1m_19580101_19581231_grid_T.nc'
df = xr.open_dataset(root_folder + fname1)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'                     

ds = xr.open_dataset(mom_dir+fname1)

# create a mask file first
# mask = np.ones(df.thetao[0,:,:,:].shape) 
# mask[np.where(df.thetao[0,:,:,:] == 0.0)] = 0
# ds2 = xr.Dataset({
#             "mask": (["deptht", "y_grid_T", "x_grid_T"], mask),
#         },
#         coords={
#             "nav_lon_grid_T": (["y_grid_T", "x_grid_T"],
#                                np.copy(df.nav_lon_grid_T)),
#             "nav_lat_grid_T": (["y_grid_T", "x_grid_T"],
#                                np.copy(df.nav_lat_grid_T)),
#             "depth":(["deptht"], np.copy(df.deptht) ),
#         },
#     )
# ds2.to_netcdf(mom_dir+'nemo0_25degree_mask.nc')

# mask = np.ones(ds.temp[0,:,:,:].shape) 
# mask[np.where(ds.temp[0,:,:,:] == 0.0)] = 0
# ds3 = xr.Dataset({
#             "mask": (["deptht", "y", "x"], mask),
#         },
#         coords={
#             "lon": (["y", "x"], np.copy(ds.lon)),
#             "lat": (["y", "x"], np.copy(ds.lat)),
#             "depth":(["deptht"], np.copy(ds.deptht) ),
#         },
#     )
# ds3.to_netcdf(mom_dir+'nemo0_25degree_mask.nc')


# ---------- define a domain target on MOM grid ---------------------
# Nx and Ny should be the hgrid sizes so double the model grid
Nx=4000
Ny=3500
domain = obc_segment('domain',path_regional_grid,istart=0,iend=Nx,jstart=0,  jend=Ny)

# ---------- define variables on each segment ------------------
north = obc_segment('segment_001',path_regional_grid,istart=Nx,iend=0, jstart=Ny,jend=Ny)
west  = obc_segment('segment_002',path_regional_grid,istart=0, iend=0, jstart=Ny,jend=0 )
south = obc_segment('segment_003',path_regional_grid,istart=0, iend=Nx,jstart=0, jend=0 )
east  = obc_segment('segment_004',path_regional_grid,istart=Nx,iend=Nx,jstart=0, jend=Ny)

# ---------- define variables on each segment ------------------
temp_south = obc_variable(south,'temp',geometry='surface',obctype='radiation',use_locstream=True)
temp_north = obc_variable(north,'temp',geometry='surface',obctype='radiation',use_locstream=True)
temp_west  = obc_variable(west, 'temp',geometry='surface',obctype='radiation',use_locstream=True)
temp_east  = obc_variable(east, 'temp',geometry='surface',obctype='radiation',use_locstream=True)

salt_south = obc_variable(south,'salt',geometry='surface',obctype='radiation',use_locstream=True)
salt_north = obc_variable(north,'salt',geometry='surface',obctype='radiation',use_locstream=True)
salt_west  = obc_variable(west, 'salt',geometry='surface',obctype='radiation',use_locstream=True)
salt_east  = obc_variable(east, 'salt',geometry='surface',obctype='radiation',use_locstream=True)

zeta_south = obc_variable(south,'zeta',geometry='line',obctype='flather',use_locstream=True)
zeta_north = obc_variable(north,'zeta',geometry='line',obctype='flather',use_locstream=True)
zeta_west  = obc_variable(west ,'zeta',geometry='line',obctype='flather',use_locstream=True)
zeta_east  = obc_variable(east ,'zeta',geometry='line',obctype='flather',use_locstream=True)

vel_south  = obc_vectvariable(south,'u','v',geometry='surface',obctype='radiation',use_locstream=True)
vel_north  = obc_vectvariable(north,'u','v',geometry='surface',obctype='radiation',use_locstream=True)
vel_west   = obc_vectvariable(west ,'u','v',geometry='surface',obctype='radiation',use_locstream=True)
vel_east   = obc_vectvariable(east ,'u','v',geometry='surface',obctype='radiation',use_locstream=True)

# ---------- run the regridding on the list of files ------------------
# for the first call to the regridding, we save all the interpolators
# (for each segment, and each type of variable), so we don't need to
# recompute the regridding weights for the N-1 following files
maskfile = mom_dir + 'nemo0_25degree_mask.nc' 
# maskfile = mom_dir + 'nemo0_25degree_mask_trim.nc' 

# compute angle 
# Approximate angles using centered differences in interior                     
lon = np.copy(df.nav_lon_grid_T)
lat = np.copy(df.nav_lat_grid_T)
angle = np.zeros((df.nav_lon_grid_T.shape[0],df.nav_lon_grid_T.shape[1]))
angle[:,1:-1] = np.arctan( (lat[:,2:]-lat[:,:-2]) /                             
                          ((lon[:,2:]-lon[:,:-2])*np.cos(np.deg2rad(lat[:,1:-1]))) )
# Approximate angles using side differences on left/right edges                 
angle[:,0] = np.arctan( (lat[:,1]-lat[:,0]) / ((lon[:,1]-lon[:,0])*np.cos(np.deg2rad(lat[:,0]))) )
angle[:,-1] = np.arctan( (lat[:,-1]-lat[:,-2]) /                                
                        ((lon[:,-1]-lon[:,-2])*np.cos(np.deg2rad(lat[:,-1]))) ) 


# output grid info                                                              
# df2 = xr.open_dataset(path_regional_grid)
# lon_rho = np.copy(df2['x'][1::2,1::2])                                          
# lat_rho = np.copy(df2['y'][1::2,1::2])                                          
# nj,ni = lon_rho.shape                                                           
# ds2 = df2['x'][1::2,1::2]                                                       
# ds2 = ds2.to_dataset(name='lon')                                                
# ds2['lat']=df2['y'][1::2,1::2]                                                  
# ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})                                  
#


first_call=False
# first_call=True
for ind in range(0,12):
    if first_call:
        interp_t2s_south = temp_south.interpolate_from(root_folder +
                                            fname1,'thetao',frame=ind,depthname='deptht',
                                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
        interp_t2s_west = temp_west.interpolate_from(root_folder +
                                            fname1,'thetao',frame=ind,depthname='deptht',
                                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
        interp_t2s_east = temp_east.interpolate_from(root_folder +
                                            fname1,'thetao',frame=ind,depthname='deptht',
                                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
        interp_t2s_north = temp_north.interpolate_from(root_folder +
                                            fname1,'thetao',frame=ind,depthname='deptht',
                                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
    else:
        temp_south.interpolate_from(root_folder + fname1,'thetao',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_south)
        temp_west.interpolate_from(root_folder + fname1,'thetao',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_west)
        temp_east.interpolate_from(root_folder + fname1,'thetao',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_east)
        temp_north.interpolate_from(root_folder + fname1,'thetao',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_north)

    salt_north.interpolate_from(root_folder + fname1,'so',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_north)
    salt_south.interpolate_from(root_folder + fname1,'so',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_south)
    salt_east.interpolate_from(root_folder + fname1,'so',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_east)
    salt_west.interpolate_from(root_folder + fname1,'so',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_west)

    zeta_north.interpolate_from(root_folder + fname1,'zos',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,autocrop=False,
                            interpolator=interp_t2s_north)
    zeta_south.interpolate_from(root_folder + fname1,'zos',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,autocrop=False,
                            interpolator=interp_t2s_south)
    zeta_east.interpolate_from(root_folder + fname1,'zos',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,autocrop=False,
                            interpolator=interp_t2s_east)
    zeta_west.interpolate_from(root_folder + fname1,'zos',frame=ind,depthname='deptht',
                            coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                            missing_value=0,autocrop=False,
                            interpolator=interp_t2s_west)

# ---------- list segments and variables to be written -------
    list_segments = [south,north,east,west]

    list_variables = [temp_south,temp_north,temp_west,temp_east,
                      salt_south,salt_north,salt_west,salt_east,
                      zeta_south,zeta_north,zeta_west,zeta_east]

    list_vectvariables = []

    #----------- time --------------------------------------------
    time = df.time_counter[ind]

    # ---------- write to file -----------------------------------
    fileout = mom_dir + fname1.replace('.nc','_obc.nc')
    write_obc_file(list_segments,list_variables,list_vectvariables,time,output=fileout)
    print(ind)
    first_call=False








# interp_t2s_south = temp_south.interpolate_from(mom_dir + fname1,'temp',frame=0,depthname='deptht',
#                                                coord_names=['lon','lat'],
#                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)

# temp_south.interpolate_from(mom_dir + fname1,'temp',frame=1,depthname='deptht',
#                             coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
#                             missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
#                             interpolator=interp_t2s_south)
