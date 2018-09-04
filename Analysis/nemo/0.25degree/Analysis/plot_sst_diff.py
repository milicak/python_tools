#!/usr/bin/env python
import numpy as np           
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import ESMF
from mpl_toolkits.basemap import Basemap                                        
import cartopy.crs as ccrs                                                      
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

plt.ion()


def enable_global(tlon,tlat,data):                                              
  """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)          
  tlon = tlon+abs(ma.max(tlon)); tlon=tlon+360                                    
  # stack grids side-by-side (in longitiudinal direction), so                   
  # any range of longitudes may be plotted on a world map.                      
  tlon = np.concatenate((tlon,tlon+360),1)                                      
  tlat = np.concatenate((tlat,tlat),1)                                          
  data = ma.concatenate((data,data),1)                                          
  tlon = tlon-360.                                                              
  return tlon, tlat, data  

maskfile = '/data/cetlod/inputs/nemo/orca025/eORCA025L75_new_subbasins.nc'
mask = xr.open_dataset(maskfile)
root_folder = '/prodigfs/fabric/NEMO_v6/ORCA025_vrac/'  
ext = '/OCE/Output/YE/'                

expid1 = 'eOR025L3P-IA-REF05-1980CLIMSTART'

#expid2 = 'eOR025L3P-IA-REF05-GMJD' 
#expid2 = 'eOR025L3P-IA-REF05' 
#expid2 = 'eOR025L3P-IA-REF05-GMJD' 
#expid2 = 'eOR025L3P-IA-REF05-LAP'
#expid2 = 'eOR025L3P-IA-REF05-noBBL'
expid2 = 'eOR025L3P-IA-REF05-GMMI' 
print expid2
                                                      
foldername1 = root_folder+expid1+ext                   
foldername2 = root_folder+expid2+ext                   
                                                      
fyear = 1990                                         
lyear = 1999                                          
                                                       
# first folder                                                       
sdate="%c%4.4d%c%c%c%c%c%c%c%c" % ('*',fyear,'0','1','0','1','*','_','T','*')
list1 = sorted(glob.glob(foldername1+'*'+sdate))        
for year in xrange(fyear+1,lyear+1):                  
   sdate="%c%4.4d%c%c%c%c%c%c%c%c" % ('*',year,'0','1','0','1','*','_','T','*')
   list1.extend(sorted(glob.glob(foldername1+'*'+sdate)))                      
                                                                             
# second folder                                                       
sdate="%c%4.4d%c%c%c%c%c%c%c%c" % ('*',fyear,'0','1','0','1','*','_','T','*')
list2 = sorted(glob.glob(foldername2+'*'+sdate))        
for year in xrange(fyear+1,lyear+1):                  
   sdate="%c%4.4d%c%c%c%c%c%c%c%c" % ('*',year,'0','1','0','1','*','_','T','*')
   list2.extend(sorted(glob.glob(foldername2+'*'+sdate)))                      
                                                                             
chunks = (134,240)                                                          
xr_chunks = {'x': chunks[-2], 'y': chunks[-1]}                              

data1 = xr.open_mfdataset(list1, chunks=xr_chunks)
data2 = xr.open_mfdataset(list2, chunks=xr_chunks)
Nx = data1.thetao.shape[-1]
Ny = data1.thetao.shape[-2]

# SST bias
sstmodel1 = np.copy(data1.thetao[:,0,:,:].mean(axis=0))
sstmodel2 = np.copy(data2.thetao[:,0,:,:].mean(axis=0))


[lon,lat,sstdiff] = enable_global(data1.nav_lon, data1.nav_lat, sstmodel2-sstmodel1)

plt.figure(figsize=(8, 6))
m=Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=90,projection='cyl')
#m = Basemap(width=12000000,height=8000000,                                     
#            resolution='l',projection='stere',\                                
#            lat_ts=40,lat_0=90,lon_0=0.)                                       
#m = Basemap(projection='stere',boundinglat=60,lon_0=0,resolution='l')          
m.drawcoastlines()                                                              
m.fillcontinents()                                                              
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])                          
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])                           
im1 = m.pcolormesh(np.transpose(lon),np.transpose(lat),np.transpose(np.ma.masked_invalid(sstdiff))
                  ,shading='flat',cmap='RdBu_r',vmin=-1,vmax=1,latlon=True)      
cb = m.colorbar(im1,"right", size="5%", pad="15%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance be 

titlename = expid2+'_SST_diff_years_'+np.str(fyear)+'_'+np.str(lyear)
plt.title(titlename)
plt.tight_layout()
savename = titlename+'.png'
plt.savefig(savename,dpi=300,bbox_inches='tight') 


sys.exit()

ax = plt.axes(projection=ccrs.PlateCarree())                                    
ax.coastlines()                                                                 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,                     
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')       
gl.ylabels_right = False                                                        
gl.xformatter = LONGITUDE_FORMATTER                                             
gl.yformatter = LATITUDE_FORMATTER

plt.pcolormesh(data1.nav_lon,data1.nav_lat,sstmodel-sstwoaann
               ,vmin=-5,vmax=5,cmap='RdBu_r',transform=ccrs.PlateCarree());plt.colorbar()

plt.colorbar(ax=ax, shrink=.75) 

plt.figure(figsize=(8, 6))
plt.pcolormesh(sstmodel-sstwoaann
               ,vmin=-5,vmax=5,cmap='RdBu_r');plt.colorbar()


