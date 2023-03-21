import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from HCtFlood.kara import flood_kara
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
from scipy.interpolate import interp1d
plt.ion()

root_folder = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/monthly/'

utss_bottom_levels= [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 
                    16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 
                    28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 
                    41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 55.0, 60.0, 65.0, 
                    70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 120.0, 140.0, 160.0,
               180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0,
               1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0 ]

utss_level_interfaces = np.roll(np.append(np.asarray(utss_bottom_levels),0),1)
interpolation_levels = 0.5*(utss_level_interfaces[1:] + utss_level_interfaces[:-1])

# if high_res true, high resolution interpolation will be used,
# If false we will use original nemo grid 
variable = 'water_level'
high_res = True
gridfile = 'utss_shyfem_esmf_meshinfo.nc'
srcgrid = ESMF.Mesh(filename=gridfile,filetype=ESMF.FileFormat.ESMFMESH)

# load nemo grid info for the T/S variables
# gr = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/nemo_BSFS/40_1d_20140716_20140812_grid_T.nc')
# gr = xr.open_dataset('/work/mi19918/Projects/uTSS/nemo_BSFS/40_1d_20140716_20140812_grid_T.nc')
gr = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/geodta/mesh_mask_bs-nrt_s5_smth_changeBosp_Marm_sill_nemo.nc')
lon_nemo = np.array(gr.nav_lon[3:30,40:85], dtype=float) 
lat_nemo = np.array(gr.nav_lat[3:30,40:85], dtype=float) 
lon_thlweg = lon_nemo.flatten()
lat_thlweg = lat_nemo.flatten()
high_res_coef1 = 1
high_res_coef2 = 1

if high_res:
    high_res_coef1 = 20
    high_res_coef2 = 27
    # interpolate nemo grid to have higher resolution
    lon_nemo = np.array(gr.nav_lon[3,40:85], dtype=float) 
    lat_nemo = np.array(gr.nav_lat[3:30,85], dtype=float) 
    Nx = lon_nemo.shape[0]
    Ny = lat_nemo.shape[0]
    lon_thlweg = lon_nemo.flatten()
    lat_thlweg = lat_nemo.flatten()
    lon = lon_thlweg
    lat = lat_thlweg
    x = np.linspace(1, lon.shape[0], num=lon.shape[0], endpoint=True)  
    y = np.linspace(1, lat.shape[0], num=lat.shape[0], endpoint=True)  
    f = interp1d(x,lon)
    g = interp1d(y,lat)
    xnew = np.linspace(1, lon.shape[0], num=high_res_coef1*lon.shape[0], endpoint=True)
    ynew = np.linspace(1, lat.shape[0], num=high_res_coef2*lat.shape[0], endpoint=True)
    lon_thlweg = f(xnew)
    lat_thlweg = g(ynew)
    lon, lat = np.meshgrid(lon_thlweg, lat_thlweg); 
    lon_thlweg = lon.flatten()
    lat_thlweg = lat.flatten()

#load the data from shyfem model
data = xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/monthly/uTSS_lobc_chunk_monthly_2020_2021.ous.nc')[variable]

coord_sys = ESMF.CoordSys.SPH_DEG
domask = True
# create locstream
locstream = ESMF.LocStream(lon_thlweg.shape[0], name="uTSS Thalweg Section", coord_sys=coord_sys)
# appoint the section locations
locstream["ESMF:Lon"] = lon_thlweg
locstream["ESMF:Lat"] = lat_thlweg
if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_thlweg.shape[0]), dtype=np.int32)



for tind in range(0,12): 
    srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
    secfield = np.zeros((Ny*high_res_coef2,Nx*high_res_coef1)) 

    # vertical level kind
    print(tind)
    srcfield.data[:] = np.copy(data[tind,:])
    srcfield.data[srcfield.data==0] = np.NaN

    # create a field on the locstream                
    dstfield = ESMF.Field(locstream, name='dstfield')
    dstfield.data[:] = 0.0                           
    
    # create an object to regrid data from the source to the destination field 
    dst_mask_values=None                                                       
    if domask:                                                                 
            dst_mask_values=np.array([0])                                      
    
    regrid = ESMF.Regrid(srcfield, dstfield,
               # filename="esmpy_example_weight_file.nc",
        # regrid_method=ESMF.RegridMethod.NEAREST_STOD,
        regrid_method=ESMF.RegridMethod.BILINEAR,
        # regrid_method=ESMF.RegridMethod.PATCH,
        unmapped_action=ESMF.UnmappedAction.IGNORE,dst_mask_values=dst_mask_values)
    
    # do the regridding from source to destination field
    dstfield = regrid(srcfield, dstfield)               
    tmp = dstfield.data
    tmp.shape = np.array([Ny*high_res_coef2,Nx*high_res_coef1])  
    secfield[:,:] = dstfield.data
 

    # create xarray dataset
    ds = xr.DataArray(secfield, coords=[lat[:,0], lon[0,:]],dims=["y", "x"]) 
    ds = ds.where(ds != 0)  
    ds = ds.to_dataset(name=variable)
    lonnew = xr.DataArray(lon, coords=[lat[:,0], lon[0,:]], dims=["y", "x"])   
    latnew = xr.DataArray(lat, coords=[lat[:,0], lon[0,:]], dims=["y", "x"])   
    ds['lon'] = lonnew 
    ds['lat'] = latnew
    outputfile1 = root_folder + '/uTSS_OBC_monthly__2020_2021_' + variable + '_nemo_grid_'  +  str(tind).zfill(2) + '.nc' 
    ds.to_netcdf(outputfile1)




