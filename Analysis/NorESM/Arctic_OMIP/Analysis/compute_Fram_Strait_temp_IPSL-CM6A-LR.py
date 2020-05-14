import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
# import ESMF
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

# if __name__ == '__main__':
#     from dask.distributed import Client, LocalCluster
#     cluster = LocalCluster(n_workers=4, memory_limit=25e9, local_dir='/cluster/work/users/cgu025/tmp_annual_mean/')
#     client = Client(cluster)

root_folder = '/tos-project1/NS9252K/CMIP6/'
omipn = 'omip1'
model = 'GFDL-CM4'
scn = 'r1i1p1f1'

fnames = root_folder+omipn+'/'+ model+'/'+scn+'/thetao_Omon_'+model+'_'+omipn+'_'+scn+'_gn*'
# fnames = root_folder+omipn+'/'+ model+'/'+scn+'/so_Omon_'+model+'_'+omipn+'_'+scn+'_gn*'

# df.thetao[0,0,310,266:279]

# 1st cycle
fnames1 = fnames+'174801-176712*'
list = sorted(glob.glob(fnames1))
CHUNKS = {'lev': 1}
# CHUNKS = {'y': 270, 'x': 360}
# df=xr.open_mfdataset(list,chunks=CHUNKS)
df=xr.open_mfdataset(list)
# GFDL indices
# Fram Strait
# i=1118:1194 j=999
# V-component
tmp = df.thetao[:,:,998,1117:1194].compute()
tmp2 = tmp.groupby('time.year').mean('time')
tmp2 = tmp2.mean('year')
# Plot figure
plt.pcolormesh(tmp.lon,-tmp.lev,tmp2,cmap='needJet2',vmin=-1.5,vmax=4,rasterized=True);plt.colorbar();plt.ylim(-4850,0);
plt.xlabel('Longitutude');
plt.ylabel('Depth [m]');
plt.savefig('paperfigs/'+model+'_'+scn+'_FramStrait_temp_1cyc.eps', bbox_inches='tight',format='eps', dpi=200)
plt.close()


# 3st cycle
fnames1 = fnames+'186801-188712*'
list = sorted(glob.glob(fnames1))
CHUNKS = {'lev': 1}
# CHUNKS = {'y': 270, 'x': 360}
# df=xr.open_mfdataset(list,chunks=CHUNKS)
df=xr.open_mfdataset(list)
# GFDL indices
# Fram Strait
# i=1118:1194 j=999
# V-component
tmp = df.thetao[:,:,998,1117:1194].compute()
tmp2 = tmp.groupby('time.year').mean('time')
tmp2 = tmp2.mean('year')
# Plot figure
plt.pcolormesh(tmp.lon,-tmp.lev,tmp2,cmap='needJet2',vmin=-1.5,vmax=4,rasterized=True);plt.colorbar();plt.ylim(-4850,0);
plt.xlabel('Longitutude');
plt.ylabel('Depth [m]');
plt.savefig('paperfigs/'+model+'_'+scn+'_FramStrait_temp_3cyc.eps', bbox_inches='tight',format='eps', dpi=200)
plt.close()

# 5st cycle
fnames1 = fnames+'198801-200712*'
list = sorted(glob.glob(fnames1))
CHUNKS = {'lev': 1}
# CHUNKS = {'y': 270, 'x': 360}
# df=xr.open_mfdataset(list,chunks=CHUNKS)
df=xr.open_mfdataset(list)
# GFDL indices
# Fram Strait
# i=1118:1194 j=999
# V-component
tmp = df.thetao[:,:,998,1117:1194].compute()
tmp2 = tmp.groupby('time.year').mean('time')
tmp2 = tmp2.mean('year')
# Plot figure
plt.pcolormesh(tmp.lon,-tmp.lev,tmp2,cmap='needJet2',vmin=-1.5,vmax=4,rasterized=True);plt.colorbar();plt.ylim(-4850,0);
plt.xlabel('Longitutude');
plt.ylabel('Depth [m]');
plt.savefig('paperfigs/'+model+'_'+scn+'_FramStrait_temp_5cyc.eps', bbox_inches='tight',format='eps', dpi=200)
plt.close()
