import numpy as np        
import cartopy.crs as ccrs 
import cartopy
import cmocean
import numpy.ma as ma     
import glob               
import xarray as xr       
import sys                
import pandas as pd       
import os                 
import matplotlib.path as mpath     
from datetime import date 

df = xr.open_dataset('/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/Exp02_0/2DArcticOcean_monthly_oceTAUY_1992_1-12.nc')    

cmap = cmocean.cm.topo 
newcmap = cmocean.tools.crop(cmap, -5000, 0, 0)  


theta = np.linspace(0, 2*np.pi, 100)                   
center, radius = [0.5, 0.5], 0.5                       
verts = np.vstack([np.sin(theta), np.cos(theta)]).T    
circle = mpath.Path(verts * radius + center)           

fig, axes = plt.subplots(figsize=(5,5))                                         
ax1 = plt.subplot2grid(shape=(1,1),loc=(0,0), colspan=1, projection=ccrs.NorthPolarStereo())                                                                                
fig.canvas.draw()                                                               
#plt.tight_layout()                                                                                                                                              
ax1.set_extent([-180, 180, 47, 90], ccrs.PlateCarree())                         
ax1.add_feature(cartopy.feature.LAND,color='grey');                             
ax1.coastlines(resolution='50m');                                               
c1 = ax1.pcolormesh(df.XC,df.YC,-ma.masked_equal(df.Depth,0),                   
                    vmin=-5000,vmax=.0,cmap=newcmap,                            
                    transform=ccrs.PlateCarree(), rasterized=True);             
#
#ax1.set_boundary(circle, transform=ax1.transAxes)   
axpos = ax1.get_position()                                                   
cbar_ax = fig.add_axes([axpos.x1+0.05,axpos.y0+0.025,0.03,axpos.height*0.97])
cbar = fig.colorbar(c1, cax=cbar_ax)                                      
plt.savefig('mitgcm_depth.png', bbox_inches='tight',format='png',dpi=300)    


