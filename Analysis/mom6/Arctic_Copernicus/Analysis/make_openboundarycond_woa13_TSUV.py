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
                                                           
momgrd = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/ocean_hgrid.nc'
woatemp = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/woa13_5564_t01_01.nc'
woasalt = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT//woa13_5564_s01_01.nc'

# ---------- define a domain target on MOM grid ---------------------
# Nx and Ny should be the hgrid sizes so double the model grid
Nx=4000
Ny=3500
# domain = obc_segment('domain',path_regional_grid,istart=0,iend=Nx,jstart=0,  jend=Ny)

# ---------- define variables on each segment ------------------
north = obc_segment('segment_001',momgrd,istart=Nx,iend=0, jstart=Ny,jend=Ny)
west  = obc_segment('segment_002',momgrd,istart=0, iend=0, jstart=Ny,jend=0 )
south = obc_segment('segment_003',momgrd,istart=0, iend=Nx,jstart=0, jend=0 )
east  = obc_segment('segment_004',momgrd,istart=Nx,iend=Nx,jstart=0, jend=Ny)

# ---------- define variables on each segment ------------------
temp_south = obc_variable(south,'temp',geometry='surface',obctype='radiation',use_locstream=True)
temp_north = obc_variable(north,'temp',geometry='surface',obctype='radiation',use_locstream=True)
temp_west  = obc_variable(west, 'temp',geometry='surface',obctype='radiation',use_locstream=True)
temp_east  = obc_variable(east, 'temp',geometry='surface',obctype='radiation',use_locstream=True)

salt_south = obc_variable(south,'salt',geometry='surface',obctype='radiation',use_locstream=True)
salt_north = obc_variable(north,'salt',geometry='surface',obctype='radiation',use_locstream=True)
salt_west  = obc_variable(west, 'salt',geometry='surface',obctype='radiation',use_locstream=True)
salt_east  = obc_variable(east, 'salt',geometry='surface',obctype='radiation',use_locstream=True)

zeta_south = obc_variable(south,'zeta',geometry='line',obctype='flather')
zeta_north = obc_variable(north,'zeta',geometry='line',obctype='flather')
zeta_west  = obc_variable(west ,'zeta',geometry='line',obctype='flather')
zeta_east  = obc_variable(east ,'zeta',geometry='line',obctype='flather')

vel_south  = obc_vectvariable(south,'u','v',geometry='surface',obctype='radiation')
vel_north  = obc_vectvariable(north,'u','v',geometry='surface',obctype='radiation')
vel_west   = obc_vectvariable(west ,'u','v',geometry='surface',obctype='radiation')
vel_east   = obc_vectvariable(east ,'u','v',geometry='surface',obctype='radiation')

days_in_month = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

# we use the depth from the temp field to create the depth vector for velocities
domain = obc_segment('domain', momgrd,istart=0,iend=Nx,jstart=0,  jend=Ny)
temp_domain = obc_variable(domain,'temp',geometry='surface',obctype='radiation')
aa=xr.open_dataset(woatemp,decode_times=False)
depth_woa=np.copy(aa.depth)

# loop over months
for kt in np.arange(12):
    # pick the correct filename
    mm=str(kt+1).zfill(2)
    woatemp = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/woa13_5564_t' + mm + '_01.nc'
    woasalt = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/woa13_5564_s' + mm + '_01.nc'
 
    # ---------- interpolate T/S from WOA file
    temp_south.interpolate_from( woatemp,'t_an',depthname='depth')
    temp_north.interpolate_from( woatemp,'t_an',depthname='depth')
    temp_west.interpolate_from(  woatemp,'t_an',depthname='depth')
    temp_east.interpolate_from(  woatemp,'t_an',depthname='depth')
    
    salt_south.interpolate_from( woasalt,'s_an',depthname='depth')
    salt_north.interpolate_from( woasalt,'s_an',depthname='depth')
    salt_west.interpolate_from(  woasalt,'s_an',depthname='depth')
    salt_east.interpolate_from(  woasalt,'s_an',depthname='depth')
    
    # ---------- set constant value for SSH and velocities -------
    zeta_south.set_constant_value(0.0)
    zeta_north.set_constant_value(0.0)
    zeta_west.set_constant_value(0.0)
    zeta_east.set_constant_value(0.0)

    vel_south.set_constant_value(0.,0.,depth_vector=depth_woa)
    vel_north.set_constant_value(0.,0.,depth_vector=depth_woa)
    vel_west.set_constant_value(0.,0.,depth_vector=depth_woa)
    vel_east.set_constant_value(0.,0.,depth_vector=depth_woa)
    
    # ---------- list segments and variables to be written -------
    list_segments = [north,south,west,east]

    list_variables = [temp_south,temp_north,temp_west,temp_east, \
                      salt_south,salt_north,salt_west,salt_east, \
                      zeta_south,zeta_north,zeta_west,zeta_east ]

    list_vectvariables = [vel_south,vel_north,vel_west,vel_east]

    #----------- time --------------------------------------------
    middle_of_current_month = days_in_month[:kt].sum() + 0.5 * days_in_month[kt]
    time = timeobject(middle_of_current_month)
    time.units = 'days since 1900-01-01'
    time.calendar = 'gregorian'

    # ---------- write to file -----------------------------------
    write_obc_file(list_segments,list_variables,list_vectvariables,time,output='/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT//obc_woa13_m' + mm + '.nc')

