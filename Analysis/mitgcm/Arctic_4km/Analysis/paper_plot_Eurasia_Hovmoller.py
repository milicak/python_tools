import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
from eofs.xarray import Eof

df0 = xr.open_dataset('Exp02_0_Eurasia_TS_mom6.nc')
df1 = xr.open_dataset('Exp02_1_Eurasia_TS_mom6.nc')
df3 = xr.open_dataset('Exp02_3_Eurasia_TS_mom6.nc')

fig, axes = plt.subplots(figsize=(6,9))
ax1 = plt.subplot2grid(shape=(4,1),loc=(0,0),rowspan=2)
ax2 = plt.subplot2grid(shape=(4,1),loc=(2,0),rowspan=2)
plt.tight_layout()

im2 = ax1.pcolormesh(df0.time,df0.depth,
               np.transpose(np.copy(df1.temp_Eurasia-df0.temp_Eurasia)),
               cmap='RdBu_r',vmin=-1,vmax=1);

ax1.text(df.time[4],-350,'a)',fontsize=14)
ax1.set_ylim(-4100,0)
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.015,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
# cbar.set_label('$^\circ$C',rotation=0,y=1.02,labelpad=-40,fontsize=14)
ax1.set_ylabel(r'Depth [m]',fontsize = 14.0)


im2 = ax2.pcolormesh(df0.time,df0.depth,
               np.transpose(np.copy(df3.temp_Eurasia-df0.temp_Eurasia)),
               cmap='RdBu_r',vmin=-1,vmax=1);

ax2.text(df.time[4],-350,'b)',fontsize=14)
ax2.set_ylim(-4100,0)
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.015,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
# cbar.set_label('$^\circ$C',rotation=0,y=1.02,labelpad=-40,fontsize=14)
ax2.set_ylabel(r'Depth [m]',fontsize = 14.0)
ax2.set_xlabel(r'year',fontsize = 14.0)


plt.savefig('paperfigs/Eurasia_temp_howmoller.png', bbox_inches='tight',format='png',dpi=300)
