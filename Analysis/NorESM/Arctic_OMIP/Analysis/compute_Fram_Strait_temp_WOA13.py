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


fnames1 = '/cluster/work/users/milicak/OMIPs_diag/woa13_decav_t00_04v2.nc'

list = sorted(glob.glob(fnames1))
CHUNKS = {'lev': 1}
# CHUNKS = {'y': 270, 'x': 360}
# df=xr.open_mfdataset(list,chunks=CHUNKS)
df = xr.open_mfdataset(list,decode_times=False)
# WOA13 indices
# Fram Strait
# i=1118:1194 j=999
# V-component
tmp2 = df.t_an[0,:,676,643:768].compute()
# Plot figure
plt.pcolormesh(df.lon[643:768],-df.depth,tmp2,cmap='needJet2',vmin=-1.5,vmax=4,rasterized=True);plt.colorbar();plt.ylim(-5550,0);
plt.xlabel('Longitutude');
plt.ylabel('Depth [m]');
plt.savefig('paperfigs/WOA13_FramStrait_temp.eps', bbox_inches='tight',format='eps', dpi=200)
plt.close()



