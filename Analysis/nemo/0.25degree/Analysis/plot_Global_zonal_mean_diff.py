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
root_folder = '/prodigfs/fabric/NEMO_v6/ORCA025_vrac/'  
ext = '/OCE/Output/YE/'                

#expid1 = 'eOR025L3P-IA-REF05' 
expid1 = 'eOR025L3P-IA-REF05-1980CLIMSTART'

#expid2 = 'eOR025L3P-IA-REF05-GMJD' 
expid2 = 'eOR025L3P-IA-REF05-GMMI' 
#expid2 = 'eOR025L3P-IA-REF05-LAP'
#expid2 = 'eOR025L3P-IA-REF05-noBBL'
print expid2
                                                      
foldername1 = root_folder+expid1+ext                   
foldername2 = root_folder+expid2+ext                   
                                                      
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

data1 = xr.open_mfdataset(list1, chunks=xr_chunks)
data2 = xr.open_mfdataset(list2, chunks=xr_chunks)
Nx = data1.thetao.shape[-1]
Ny = data1.thetao.shape[-2]
# mask variable
# Global 
maskvarG = mask.glomsk.where(mask.glomsk==1)
# Atlantic
maskvarA = mask.atlmsk.where(mask.atlmsk==1)
# Pacific
maskvarP = mask.pacmsk.where(mask.pacmsk==1)
# Indian
maskvarI = mask.indmsk.where(mask.indmsk==1)

# Global zonal mean bias
data1 = data1*maskvarG
data2 = data2*maskvarG
lat = np.copy(grid.nav_lat*maskvarG)
#lat[lat==0] = np.nan
lat = np.nanmean(lat,axis=1)
lat[np.isnan(lat)] = 0

datznl1 = np.copy(data1.thetao.mean(dim=['time_counter','x']));
datznl2 = np.copy(data2.thetao.mean(dim=['time_counter','x']));
#sys.exit()

cmap1 = cm.get_cmap("RdBu_r")
#cmap1 = cm.get_cmap("RdBu_r",lut=12) # how many layers 
cmap1.set_bad("grey")

#plt.figure()
#plt.pcolormesh(np.linspace(1,Ny,Ny),-data1.deptht,(datznl-woaznl)
#               ,vmin=-2,vmax=2,cmap='jet');plt.colorbar()

plt.figure()
plt.pcolormesh(lat,-data1.deptht,(datznl2-datznl1)
               ,vmin=-1,vmax=1,cmap=cmap1);plt.colorbar()

titlename = expid2+'_Global_diff_years_'+np.str(fyear)+'_'+np.str(lyear)
plt.title(titlename)
plt.tight_layout()
savename = titlename+'.png'
plt.savefig(savename,dpi=300,bbox_inches='tight') 
