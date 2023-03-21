'''This is code is written using
   https://github.com/ESMG/PyCNAL_regridding/blob/master/examples/SODA3.3.1/Creating_Initial_and_Boundary_conditions_from_SODA3.py'''
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
ds = xr.open_dataset(path_regional_grid)
ls1 = sorted(glob.glob(root_folder+'*grid_T*'))
xstr = 305
# xend = 314
xend = 366

# ds = xr.open_dataset(mom_dir+fname1)
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
Ny=2600
# domain = obc_segment('domain',path_regional_grid,istart=0,iend=Nx,jstart=0,  jend=Ny)

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


# first_call=False
first_call=True
days = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

mid_days = np.zeros(days.shape)
mid_days[0] = 0.5*days[0]
for ind in range(1,12):
    mid_days[ind] = 0.5*days[ind] + days[:ind].sum()



for ind0 in range(xstr,xend):
    fname1 = ls1[ind0][-44:]
    fname2 = fname1[:-3]+'_rotated.nc'
    for ind in range(0,12):
        df = xr.open_dataset(root_folder + fname1)
        tail = '_obc_' + str(ind).zfill(2) + '.nc'
        fileout = mom_dir + fname1.replace('.nc',tail)
        if not os.path.isfile(fileout):
            if first_call:
                interp_t2s_south = temp_south.interpolate_from(mom_dir +
                                                fname1,'thetao',frame=ind,depthname='deptht',
                                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
                interp_t2s_west = temp_west.interpolate_from(mom_dir +
                                                fname1,'thetao',frame=ind,depthname='deptht',
                                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
                interp_t2s_east = temp_east.interpolate_from(mom_dir +
                                                fname1,'thetao',frame=ind,depthname='deptht',
                                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
                interp_t2s_north = temp_north.interpolate_from(mom_dir +
                                                fname1,'thetao',frame=ind,depthname='deptht',
                                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
            else:
                temp_south.interpolate_from(mom_dir + fname1,'thetao',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_south)
                temp_west.interpolate_from(mom_dir + fname1,'thetao',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_west)
                temp_east.interpolate_from(mom_dir + fname1,'thetao',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_east)
                temp_north.interpolate_from(mom_dir + fname1,'thetao',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_north)
    
            salt_north.interpolate_from(mom_dir + fname1,'so',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_north)
            salt_south.interpolate_from(mom_dir + fname1,'so',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_south)
            salt_east.interpolate_from(mom_dir + fname1,'so',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_east)
            salt_west.interpolate_from(mom_dir + fname1,'so',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator=interp_t2s_west)
    
            zeta_north.interpolate_from(mom_dir + fname1,'zos',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,autocrop=False,
                                interpolator=interp_t2s_north)
            zeta_south.interpolate_from(mom_dir + fname1,'zos',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,autocrop=False,
                                interpolator=interp_t2s_south)
            zeta_east.interpolate_from(mom_dir + fname1,'zos',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,autocrop=False,
                                interpolator=interp_t2s_east)
            zeta_west.interpolate_from(mom_dir + fname1,'zos',frame=ind,depthname='deptht',
                                coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
                                missing_value=0,autocrop=False,
                                interpolator=interp_t2s_west)
    
            if first_call:
                interp_u2s_south, interp_v2s_south = vel_south.interpolate_from(mom_dir +
                                                fname2,'ue','vn',frame=ind,depthname='deptht',
                                                coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                                coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
                interp_u2s_west, interp_v2s_west   = vel_west.interpolate_from(mom_dir +
                                                fname2,'ue','vn',frame=ind,depthname='deptht',
                                                coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                                coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
                interp_u2s_east, interp_v2s_east   = vel_east.interpolate_from(mom_dir +
                                                fname2,'ue','vn',frame=ind,depthname='deptht',
                                                coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                                coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
                interp_u2s_north, interp_v2s_north = vel_north.interpolate_from(mom_dir +
                                                fname2,'ue','vn',frame=ind,depthname='deptht',
                                                coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                                coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)  
            else:
                vel_north.interpolate_from(mom_dir + fname2,'ue','vn',frame=ind,depthname='deptht',
                                    coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                    coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                    missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                    interpolator_u=interp_u2s_north,interpolator_v=interp_v2s_north)
                vel_south.interpolate_from(mom_dir + fname2,'ue','vn',frame=ind,depthname='deptht',
                                    coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                    coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                    missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                    interpolator_u=interp_u2s_south,interpolator_v=interp_v2s_south)
                vel_east.interpolate_from(mom_dir + fname2,'ue','vn',frame=ind,depthname='deptht',
                                    coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                    coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                    missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                    interpolator_u=interp_u2s_east,interpolator_v=interp_v2s_east)
                vel_west.interpolate_from(mom_dir + fname2,'ue','vn',frame=ind,depthname='deptht',
                                    coord_names_u=['nav_lon_grid_T','nav_lat_grid_T'],
                                    coord_names_v=['nav_lon_grid_T','nav_lat_grid_T'],
                                    missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
                                    interpolator_u=interp_u2s_west,interpolator_v=interp_v2s_west)


            # DONT FORGET PYCNAL regridding already performed apply rotation transpose

    # ---------- list segments and variables to be written -------
            list_segments = [south,north,east,west]
    
            list_variables = [temp_south,temp_north,temp_west,temp_east,
                              salt_south,salt_north,salt_west,salt_east,
                              zeta_south,zeta_north,zeta_west,zeta_east]
    
           # list_vectvariables = []
            list_vectvariables = [vel_south,vel_north,vel_west,vel_east]
    
            #----------- time --------------------------------------------
            # time = df.time_counter[ind]
            time = timeobject(mid_days[ind])
            time.units = 'days since 1958-01-01'
            time.calendar = 'NoLeap'
            # time.attrs['units']='days since 1900-01-01'
            # time.attrs['calendar'] = 'gregorian'
    
            # ---------- write to file -----------------------------------
            write_obc_file(list_segments,list_variables,list_vectvariables,time,output=fileout)
            print(fname1)
            print(ind)
            first_call=False








# interp_t2s_south = temp_south.interpolate_from(mom_dir + fname1,'temp',frame=0,depthname='deptht',
#                                                coord_names=['lon','lat'],
#                                                missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False)

# temp_south.interpolate_from(mom_dir + fname1,'temp',frame=1,depthname='deptht',
#                             coord_names=['nav_lon_grid_T','nav_lat_grid_T'],
#                             missing_value=0,maskfile=maskfile,maskvar='mask',autocrop=False,
#                             interpolator=interp_t2s_south)
