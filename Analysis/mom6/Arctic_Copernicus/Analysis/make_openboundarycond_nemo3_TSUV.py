'''This is code is written using
   https://github.com/ESMG/PyCNAL_regridding/blob/master/examples/SODA3.3.1/Creating_Initial_and_Boundary_conditions_from_SODA3.py'''
import os                                                  
import numpy as np                                         
import xesmf as xe                                         
import xesmf
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


def open_grid(path,decode_times=False):
    """Return a grid object containing staggered grid locations"""
    grid={}
    grid['ds']=xr.open_dataset(path,decode_times=False)
    grid['ds']=grid['ds'].drop_dims(['ny','nx'])
    grid['ds']=grid['ds'].drop_vars(['tile'])
    grid['nyp']=grid['ds'].nyp.data[-1]+1
    grid['nxp']=grid['ds'].nxp.data[-1]+1
    nxp=grid['nxp'];nyp=grid['nyp']
    grid['h'] = grid['ds'].isel(nxp=slice(1,nxp+1,2),nyp=slice(1,nyp+1,2))
    #The q grid is not symmetric, but Cu and Cv are
    grid['q'] = grid['ds'].isel(nxp=slice(2,nxp+1,2),nyp=slice(2,nyp+1,2))
    grid['Cu'] = grid['ds'].isel(nxp=slice(0,nxp+1,2),nyp=slice(1,nyp+1,2))
    grid['Cv'] = grid['ds'].isel(nxp=slice(1,nxp+1,2),nyp=slice(0,nyp+1,2))
    return grid

                                                           
root_folder = '/data/inputs/metocean/historical/model/ocean/MERCATOR/CMEMS/reanalysis/day/'
fname1 = '1993/01/grepv2_daily_mnstd_19930101.nc'
df = xr.open_dataset(root_folder + fname1)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'                     
regional_grid=open_grid(path_regional_grid)


# Define boundaries for regional domain model
ds_regional=regional_grid['ds']
# northern boundary
north = xr.Dataset()
north['lon'] = ds_regional['x'].isel(nyp=-1)
north['lat'] = ds_regional['y'].isel(nyp=-1)
# southern boundary
south = xr.Dataset()
south['lon'] = ds_regional['x'].isel(nyp=0)
south['lat'] = ds_regional['y'].isel(nyp=0)
# western boundary
west = xr.Dataset()
west['lon'] = ds_regional['x'].isel(nxp=0)
west['lat'] = ds_regional['y'].isel(nxp=0)
# eastern boundary
east = xr.Dataset()
east['lon'] = ds_regional['x'].isel(nxp=-1)
east['lat'] = ds_regional['y'].isel(nxp=-1)


ls1 = sorted(glob.glob('/data/inputs/metocean/historical/model/ocean/MERCATOR/CMEMS/reanalysis/day/1993/*/grepv2_daily_mnstd*'))
dfm = xr.open_mfdataset(ls1)
dfm = dfm.resample(time='5D').mean('time')
dfm = dfm.rename({'longitude':'lon','latitude':'lat'})

ds_uvq_r = xr.Dataset({'u':uq_rot,'v':vq_rot},coords={'time':ds_uvq.time,'lon':ds_uvq.lon,'lat':ds_uvq.lat})


# Calculate remapping weights
# Using nearest neighbor - other options could be used here , e.g. bilinear.
regrid_north_uv = xesmf.Regridder(dfm, north, 'nearest_s2d', 
                               locstream_out=True, periodic=False, filename='regrid_north_uv.nc',reuse_weights=True)


#Note that parent grid uv values are symmetric
path_parent_grid='/net2/mjh/ipynb/OM4_025/c192_mosaic/ocean_hgrid.nc'
parent_grid=open_grid(path_parent_grid)


ls1 = sorted(glob.glob('/data/inputs/metocean/historical/model/ocean/MERCATOR/CMEMS/reanalysis/day/*/*/grepv2_daily_mnstd*'))
xstr = 0
xend = 9861

# create a mask file first
# mask = np.ones(df.thetao_mean[0,:,:,:].shape) 
# mask[np.where(np.isnan(df.thetao_mean[0,:,:,:]))]=0
# ds2 = xr.Dataset({
#             "mask": (["depth", "latitude", "longitude"], mask),
#         },
#         coords={
#             "longitude": (["longitude"],
#                                np.copy(df.longitude)),
#             "latitude": (["latitude"],
#                                np.copy(df.latitude)),
#             "depth":(["depth"], np.copy(df.depth) ),
#         },
#     )
# ds2.to_netcdf(mom_dir+'nemo_cmems_reanalysis_mask.nc')

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
maskfile = mom_dir + 'nemo_cmems_reanalysis_mask.nc' 

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

# set variables
depth_var = 'depth'
lon_var = 'longitude'
lat_var = 'latitude'
theta_var = 'thetao_mean'
so_var = 'so_mean'
zos_var = 'zos_mean'
ue_var = 'uo_mean' 
vn_var = 'vo_mean' 
tail = '_obc.nc' 

for ind0 in range(xstr,xend):
    fname1 = ls1[ind0]
    fname2 = ls1[ind0]
    fileout = mom_dir + fname1[-30:].replace('.nc',tail)
    if not os.path.isfile(fileout):
        if first_call:
            interp_t2s_south = temp_south.interpolate_from(
                                            fname1,theta_var,frame=0,depthname=depth_var,
                                            coord_names=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
            interp_t2s_west = temp_west.interpolate_from(
                                            fname1,theta_var,frame=0,depthname=depth_var,
                                            coord_names=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
            interp_t2s_east = temp_east.interpolate_from(
                                            fname1,theta_var,frame=0,depthname=depth_var,
                                            coord_names=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
            interp_t2s_north = temp_north.interpolate_from(
                                            fname1,theta_var,frame=0,depthname=depth_var,
                                            coord_names=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
        else:
            temp_south.interpolate_from(
                                        fname1,theta_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_south)
            temp_west.interpolate_from(
                                       fname1,theta_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_west)
            temp_east.interpolate_from(
                                       fname1,theta_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_east)
            temp_north.interpolate_from(
                                        fname1,theta_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_north)
    
        salt_north.interpolate_from(
                                    fname1,so_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_north)
        salt_south.interpolate_from(
                                    fname1,so_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_south)
        salt_east.interpolate_from(
                                   fname1,so_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_east)
        salt_west.interpolate_from(
                                   fname1,so_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                            interpolator=interp_t2s_west)
    
        zeta_north.interpolate_from(
                                    fname1,zos_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,autocrop=False,
                            interpolator=interp_t2s_north)
        zeta_south.interpolate_from(
                                    fname1,zos_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,autocrop=False,
                            interpolator=interp_t2s_south)
        zeta_east.interpolate_from(
                                   fname1,zos_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,autocrop=False,
                            interpolator=interp_t2s_east)
        zeta_west.interpolate_from(
                                   fname1,zos_var,frame=0,depthname=depth_var,
                            coord_names=[lon_var,lat_var],
                            missing_value=None,autocrop=False,
                            interpolator=interp_t2s_west)
    
        if first_call:
            interp_u2s_south, interp_v2s_south = vel_south.interpolate_from(
                                            fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                            coord_names_u=[lon_var,lat_var],
                                            coord_names_v=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
            interp_u2s_west, interp_v2s_west   = vel_west.interpolate_from(
                                            fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                            coord_names_u=[lon_var,lat_var],
                                            coord_names_v=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
            interp_u2s_east, interp_v2s_east   = vel_east.interpolate_from(
                                            fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                            coord_names_u=[lon_var,lat_var],
                                            coord_names_v=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
            interp_u2s_north, interp_v2s_north = vel_north.interpolate_from(
                                            fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                            coord_names_u=[lon_var,lat_var],
                                            coord_names_v=[lon_var,lat_var],
                                            missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False)  
        else:
            vel_north.interpolate_from(
                                       fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                coord_names_u=[lon_var,lat_var],
                                coord_names_v=[lon_var,lat_var],
                                missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator_u=interp_u2s_north,interpolator_v=interp_v2s_north)
            vel_south.interpolate_from(
                                       fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                coord_names_u=[lon_var,lat_var],
                                coord_names_v=[lon_var,lat_var],
                                missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator_u=interp_u2s_south,interpolator_v=interp_v2s_south)
            vel_east.interpolate_from(
                                      fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                coord_names_u=[lon_var,lat_var],
                                coord_names_v=[lon_var,lat_var],
                                missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator_u=interp_u2s_east,interpolator_v=interp_v2s_east)
            vel_west.interpolate_from(
                                      fname2,ue_var,vn_var,frame=0,depthname=depth_var,
                                coord_names_u=[lon_var,lat_var],
                                coord_names_v=[lon_var,lat_var],
                                missing_value=None,maskfile=maskfile,maskvar='mask',autocrop=False,
                                interpolator_u=interp_u2s_west,interpolator_v=interp_v2s_west)


        # DONT FORGET PYCNAL regridding already performed apply rotation transpose
    # -------- list segments and variables to be written -------
        list_segments = [south,north,east,west]
    
        list_variables = [temp_south,temp_north,temp_west,temp_east,
                          salt_south,salt_north,salt_west,salt_east,
                          zeta_south,zeta_north,zeta_west,zeta_east]
    
       # list_vectvariables = []
        list_vectvariables = [vel_south,vel_north,vel_west,vel_east]
    
        #----------- time --------------------------------------------
        # time = df.time_counter[ind]
        time = timeobject(ind0)
        time.units = 'days since 1993-01-01'
        time.calendar = 'gregorian'
        # time.calendar = 'NoLeap'
        # time.attrs['units']='days since 1900-01-01'
        # time.attrs['calendar'] = 'gregorian'
    
        # ---------- write to file -----------------------------------
        write_obc_file(list_segments,list_variables,list_vectvariables,time,output=fileout)
        print(fname1)
        first_call=False








