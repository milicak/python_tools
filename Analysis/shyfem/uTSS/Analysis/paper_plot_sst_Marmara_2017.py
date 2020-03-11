import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
# import ESMF
from mpl_toolkits.basemap import Basemap                                            
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()

def plot_fullmap(df,variable):
    var = df[variable][0,:,0]
    plt.figure()
    m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc',\
                lat_0=40.,lon_0=20.,lat_ts=20.)
    
    m.drawcoastlines(linewidth=0.2)
    m.fillcontinents(color='grey')
    m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0])
    m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])
    longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))
    im1=plt.tripcolor(longitude,latitude,df.element_index,
                    np.ma.masked_equal(var,0),cmap='needJet2',vmin=22,vmax=28,shading='gouraud')
    cb = m.colorbar(im1,"right", size="5%", pad="10%")
    savename = 'paperfigs/uTSS_SST_Marmara_'+\
                np.str(np.copy(df.time.dt.day))[1:3]+'_'+\
                np.str(np.copy(df.time.dt.month))[1:2].zfill(2)+'_'+\
                np.str(np.copy(df.time.dt.year))[1:5]+'.png'
    plt.savefig(savename, bbox_inches='tight',format='png', dpi=300)
    plt.close()

def plot_zoommap(df,variable):
    var = df[variable][0,:,0]
    plt.figure()
    m = Basemap(llcrnrlon=26.0,llcrnrlat=40.0,urcrnrlon=30.,urcrnrlat=41.2,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc',\
                lat_0=40.,lon_0=20.,lat_ts=20.)
    
    m.drawcoastlines(linewidth=0.2)
    m.fillcontinents(color='grey')
    m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0])
    m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])
    longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))
    im1=plt.tripcolor(longitude,latitude,df.element_index,
                    np.ma.masked_equal(var,0),cmap='needJet2',vmin=22,vmax=28,shading='gouraud')
    cb = m.colorbar(im1,"right", size="5%", pad="10%")
    savename = 'paperfigs/uTSS_SST_Marmara_'+\
                np.str(np.copy(df.time.dt.day))[1:3]+'_'+\
                np.str(np.copy(df.time.dt.month))[1:2].zfill(2)+'_'+\
                np.str(np.copy(df.time.dt.year))[1:5]+'_zoom.png'
    plt.savefig(savename, bbox_inches='tight',format='png', dpi=300)
    plt.close()


root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

# expid = 'Exp01.2'
# expid = 'Exp_20160101'
expid = 'Exp_2016_analysis_newTSIC'

fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_0590.nos.nc'

df = xr.open_dataset(fname)
df['element_index'] -= 1

plot_zoommap(df,'temperature')


