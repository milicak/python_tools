
import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import rc
from needJet2 import shfn
import ESMF
from mpl_toolkits.basemap import Basemap                                            
#import cartopy.crs as ccrs                                                          
import scipy.io as sio
import geopy.distance


aa = sio.loadmat('thalweg_xy.mat')

secdist = np.zeros(aa['lonthalweg'].shape[0])
for kind in range(0,aa['lonthalweg'].shape[0]-1):
    coords_1 = (np.copy(aa['latthalweg'][kind]),np.copy(aa['lonthalweg'][kind]))
    coords_2 = (np.copy(aa['latthalweg'][kind+1]),np.copy(aa['lonthalweg'][kind+1]))
    secdist[kind+1]=geopy.distance.distance(coords_1, coords_2).m


root_folder  = '/work/mi19918/Projects/'                                         
project_name = 'uTSS'                                                            
                                                                                 
expid = 'Exp_20160101'                                                           
#expid = 'Exp01.2'                                                               
                                                                                 
#fname = root_folder+project_name+'/'+expid+'/'+'uTSSm0_ous.nc'                  
fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_0000.nos.nc' 
                                                                                 
data = xr.open_dataset(fname)                                                    
gridfile = 'utss_shyfem_esmf_meshinfo.nc'

coord_sys=ESMF.CoordSys.SPH_DEG                                                 
domask=True                                                                     
# create locstream                                                              
locstream = ESMF.LocStream(aa['lonthalweg'].shape[0], name="TSS Thalweg Section", coord_sys=coord_sys)
# appoint the section locations                                                 
locstream["ESMF:Lon"] = aa['lonthalweg'][:,0]                                               
locstream["ESMF:Lat"] = aa['latthalweg'][:,0]                                       
if domask:                                                                      
    locstream["ESMF:Mask"] = np.array(np.ones(aa['lonthalweg'].shape[0]), dtype=np.int32)


# Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
srcgrid = ESMF.Mesh(filename=gridfile,filetype=ESMF.FileFormat.ESMFMESH)                         
srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER) 
dstfield = ESMF.Field(locstream, name='dstfield')  
# create an object to regrid data from the source to the destination field 
dst_mask_values=None                                                       
if domask:                                                                 
    dst_mask_values=np.array([0])                                      

secfield = np.zeros((aa['lonthalweg'].shape[0],data.salinity.shape[2])) 


for kind in range(0,data.salinity.shape[2]):
    print 'indice = ', kind 
    srcfield.data[:] = data.salinity[0,:,kind]     
                                                   
    # create a field on the locstream                  
    dstfield.data[:] = 0.0                             
                                                                               
    regrid = ESMF.Regrid(srcfield, dstfield,                                   
                        # regrid_method=ESMF.RegridMethod.BILINEAR,              
                        regrid_method=ESMF.RegridMethod.NEAREST_STOD,              
                        unmapped_action=ESMF.UnmappedAction.IGNORE,            
                        dst_mask_values=dst_mask_values)                       
    
    # do the regridding from source to destination field   
    dstfield = regrid(srcfield, dstfield)                  
    secfield[:,kind] = dstfield.data                       
    


plt.pcolormesh(np.cumsum(secdist),-np.cumsum(data.level,axis=0),ma.masked_equal(np.transpose(secfield),0),cmap='jet');

