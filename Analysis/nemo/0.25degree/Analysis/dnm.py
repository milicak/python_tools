#!/usr/bin/env python
import numpy as np           
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ESMF

plt.ion()

maskfile = '/data/cetlod/inputs/nemo/orca025/eORCA025L75_new_subbasins.nc'
mask = xr.open_dataset(maskfile)
grid = xr.open_dataset('/prodigfs/fabric/NEMO_v6/mesh_masks/mesh_mask_eORCA025.nc')
root_folder1 = '/prodigfs/fabric/NEMO_v6/ORCA025_vrac/'  

#root_folder2 = '/prodigfs/fabric/NEMO_v6/ORCA025_vrac/'  
root_folder2 = '/ccc/store/cont003/thredds/personr/IGCM_OUT/ORCA025_LIM3/DEVT/ORCA025ia/'
ext1 = '/OCE/Output/YE/'                
#ext2 = '/OCE/Output/YE/'    
ext2 = '/OCE/Output/MO/'    

#expid1 = 'eOR025L3P-IA-REF05' 
expid1 = 'eOR025L3P-IA-REF05-1980CLIMSTART'

#expid2 = 'eOR025L3P-IA-REF05-GMJD' 
#expid2 = 'eOR025L3P-IA-REF05-GMMI' 
#expid2 = 'eOR025L3P-IA-REF05-LAP'
#expid2 = 'eOR025L3P-IA-REF05-noBBL'
#expid2 = 'eOR025L3P-IA-REF05-GM'
#expid2 = 'eOR025L3-IA-REF05-GMv2'
#expid2 = 'eOR025L3-IA-REF05-GMMD'
expid2 = 'eOR025L3-IA-REF05-FK'
#expid2 = 'eOR025L3-IA-REF05-GMMDv2'
print expid2
                                                      
foldername1 = root_folder1+expid1+ext1 
foldername2 = root_folder2+expid2+ext2                  
                                                      
fyear = 1990                                         
lyear = 1999                                          
                                                       
# first folder                                                       
sdate="%c%4.4d%c%c%c%c%c%c%c%c" % ('*',fyear,'0','1','0','1','*','_','T','*')
list1 = sorted(glob.glob(foldername1+'*'+sdate))        
list2 = sorted(glob.glob(foldername2+'*'+sdate))        
for year in xrange(fyear+1,lyear+1):                  
   sdate="%c%4.4d%c%c%c%c%c%c%c%c" % ('*',year,'0','1','0','1','*','_','T','*')
   list1.extend(sorted(glob.glob(foldername1+'*'+sdate)))                      
   list2.extend(sorted(glob.glob(foldername2+'*'+sdate)))                      
                                                                             

chunks = (134,240)                                                          
xr_chunks = {'x': chunks[-2], 'y': chunks[-1]}                              
xr_chunks2 = {'time_counter': 5, 'x': chunks[-2], 'y': chunks[-1]}                              

data1 = xr.open_mfdataset(list1, chunks=xr_chunks,decode_times=False)
