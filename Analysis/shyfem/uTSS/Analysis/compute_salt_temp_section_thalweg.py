import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
from mpl_toolkits.basemap import Basemap                                            
from scipy.interpolate import interp1d
import geopy.distance
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()

root_folder = '/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/'

fname = root_folder + 'uTSS_TS_monthly_clim.nc'

df = xr.open_dataset(fname)


# aa = np.loadtxt('thalwegCoords_sortednew.txt')
aa = np.loadtxt('newthalweg.txt')
gridfile = 'utss_shyfem_esmf_meshinfo.nc'
srcgrid = ESMF.Mesh(filename=gridfile,filetype=ESMF.FileFormat.ESMFMESH)

lon_thlweg = aa[:,0]
lat_thlweg = aa[:,1]
lon = lon_thlweg
lat = lat_thlweg
x = np.linspace(1, lon.shape[0], num=lon.shape[0], endpoint=True)  
f = interp1d(x,lon)
g = interp1d(x,lat)
xnew = np.linspace(1, lon.shape[0], num=2*lon.shape[0], endpoint=True)
lon_thlweg = f(xnew)
lat_thlweg = g(xnew)
dnm = np.vstack([lat_thlweg,lon_thlweg])
dist = np.zeros(lon_thlweg.shape)
for ind in range(0,lon_thlweg.shape[0]-1):
    dist[ind+1] = dist[ind]+geopy.distance.distance(dnm[:,ind+1],dnm[:,ind]).km

ds = df.mean('month')

coord_sys = ESMF.CoordSys.SPH_DEG
domask = True
# create locstream
locstream = ESMF.LocStream(lon_thlweg.shape[0], name="uTSS Thalweg Section", coord_sys=coord_sys)
# appoint the section locations
locstream["ESMF:Lon"] = lon_thlweg
locstream["ESMF:Lat"] = lat_thlweg
if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_thlweg.shape[0]), dtype=np.int32)


secfield = np.zeros((lon_thlweg.shape[0],np.copy(df.level.shape[0])))
srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)

# vertical level kind
for kind in range(0,93):
    print(kind)
    srcfield.data[:] = np.copy(ds.salinity[:,kind])
    srcfield.data[srcfield.data==0] = np.NaN

    # create a field on the locstream                
    dstfield = ESMF.Field(locstream, name='dstfield')
    dstfield.data[:] = 0.0                           
    
    # create an object to regrid data from the source to the destination field 
    dst_mask_values=None                                                       
    if domask:                                                                 
            dst_mask_values=np.array([0])                                      
    
    regrid = ESMF.Regrid(srcfield, dstfield,
        # regrid_method=ESMF.RegridMethod.NEAREST_STOD,
        regrid_method=ESMF.RegridMethod.BILINEAR,
        # regrid_method=ESMF.RegridMethod.PATCH,
        unmapped_action=ESMF.UnmappedAction.IGNORE,dst_mask_values=dst_mask_values)
    
    # do the regridding from source to destination field
    dstfield = regrid(srcfield, dstfield)               
    secfield[:,kind] = dstfield.data
 

salt_sec = xr.DataArray(secfield, dims=("dist", "zr"))

secfield = np.zeros((lon_thlweg.shape[0],np.copy(df.level.shape[0])))
srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
# vertical level kind
for kind in range(0,93):
    print(kind)
    srcfield.data[:] = np.copy(ds.temperature[:,kind])
    srcfield.data[srcfield.data==0] = np.NaN

    # create a field on the locstream                
    dstfield = ESMF.Field(locstream, name='dstfield')
    dstfield.data[:] = 0.0                           
    
    # create an object to regrid data from the source to the destination field 
    dst_mask_values=None                                                       
    if domask:                                                                 
            dst_mask_values=np.array([0])                                      
    
    regrid = ESMF.Regrid(srcfield, dstfield,
        regrid_method=ESMF.RegridMethod.BILINEAR,
        unmapped_action=ESMF.UnmappedAction.IGNORE,dst_mask_values=dst_mask_values)
    
    # do the regridding from source to destination field
    dstfield = regrid(srcfield, dstfield)               
    secfield[:,kind] = dstfield.data
 

temp_sec = xr.DataArray(secfield, dims=("dist", "zr"))

df1 = xr.Dataset({"salt_thalweg": salt_sec, "temp_thalweg": temp_sec,
                  "dist": dist , "zr": np.copy(df.level)})  

df1.to_netcdf('temp_salt_thalweg_sections.nc')

plt.figure()
plt.pcolormesh(dist,
               -df.level,np.transpose(ma.masked_invalid(secfield))
              ,cmap='needJet2');
plt.colorbar()



