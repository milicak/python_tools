import numpy as np
import xarray as xr


df1 = xr.open_dataset('Exp01_0_zonal_mean_stress.nc')
df2 = xr.open_dataset('Exp02_0_zonal_mean_stress.nc')
df3 = xr.open_dataset('Exp03_0_zonal_mean_stress.nc')

root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp01_0'
grname = root_folder+project_name+'/grid.nc'
gr = xr.open_dataset(grname)


plt.plot(df1.oceTAUX,gr.YC,'b',label='CTRL')
plt.plot(df3.oceTAUX,gr.YC,'r',label='THERMAL')
plt.legend(fontsize=14);
plt.grid()
plt.xticks(fontsize=14);
plt.yticks(fontsize=14);
plt.xlabel('Nm$^{-2}$',fontsize=14);
plt.ylabel('Lat',fontsize=14);

plt.savefig('paperfigs/zonal_wind_stress.png', bbox_inches='tight',format='png',dpi=300)
