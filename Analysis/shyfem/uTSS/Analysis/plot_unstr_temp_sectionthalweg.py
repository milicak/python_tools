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
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

# expid = 'Exp01.2'
# expid = 'Exp_20160101'
expid = 'Exp_2016_analysis'
# expid = 'Exp_2016_analysis_keps'

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

fyear = 0;
lyear = 0;
# lyear = 396;

sdate = "%c%4.4d%c" % ('*',fyear,'*')
fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
list = sorted(glob.glob(fname))
str1 = ''.join(list)
grd = xr.open_dataset(str1)
for year in xrange(fyear+1,lyear+1):
    sdate = "%c%4.4d%c" % ('*',year,'*')
    list.extend(sorted(glob.glob(root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc')))



# data = xr.open_mfdataset(list, chunks={'time':5, 'node':20000})['salinity']
data = xr.open_mfdataset(list, chunks={'time':5, 'node':20000})['temperature']
# data['element_index'] -= 1
grd['element_index'] -= 1
ds = data.mean('time')

coord_sys = ESMF.CoordSys.SPH_DEG
domask = True
# create locstream
locstream = ESMF.LocStream(lon_thlweg.shape[0], name="uTSS Thalweg Section", coord_sys=coord_sys)
# appoint the section locations
locstream["ESMF:Lon"] = lon_thlweg
locstream["ESMF:Lat"] = lat_thlweg
if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_thlweg.shape[0]), dtype=np.int32)


secfield = np.zeros((lon_thlweg.shape[0],np.copy(data.level.shape[0])))

srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)

# vertical level kind
for kind in range(0,93):
    print(kind)
    # srcfield.data[:] = np.copy(data[0,:,kind])
    srcfield.data[:] = np.copy(ds[:,kind])
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
 

plt.figure()
# triang=mtri.Triangulation(grd.longitude,grd.latitude,grd.element_index)
# plt.pcolormesh(np.linspace(1,295,295),-data.level,np.transpose(ma.masked_equal(secfield-secfieldold,0)));
plt.pcolormesh(np.linspace(1,lon_thlweg.shape[0],lon_thlweg.shape[0]),
               # -data.level,np.transpose(ma.masked_equal(secfield,0))
               -data.level,np.transpose(ma.masked_invalid(secfield))
              ,cmap='needJet2');
plt.colorbar()



