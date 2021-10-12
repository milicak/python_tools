import numpy as np
from datetime import datetime
from datetime import timedelta  
import sys
from scipy.spatial import cKDTree                                                                                                                               

startTime = datetime.now()

### to be set

date_init 	= datetime(2018,1,1,0,0,0)	
date_init_str 	= date_init.strftime('%Y%m%d %H%M%S\n')


df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/uTSS_UV_monthly_clim.nc')  

ssh	= df['water_level']
nkn     = df.longitude.shape[1]
lon_utss = np.copy(df.longitude[-1,:])
lat_utss = np.copy(df.latitude[-1,:])

# 1509 tsunami values
ds = pd.read_excel('Istanbul_tsunami.xlsx')   
lon_t = ds['Longitude']
lat_t = ds['Latitude']
ssh_t = ds['DZ 1509']
lon_t = np.array(lon_t)
lat_t = np.array(lat_t)
ssh_t = np.array(ssh_t)
aa = np.vstack((lon_t,lat_t))  
aa = np.transpose(aa)

# interpolate
ssh_utss = griddata(aa, ssh_t, (lon_utss, lat_utss), method='linear') 
ssh_utss[np.isnan(ssh_utss)==1] = 0

# select first time value
ssh	= np.copy(ssh[0,:])

# ssh = ssh + ssh_utss
ssh = ssh_utss


###########################
### WRITE TO FILE
###########################

print('write water level file..')
header = '0 2 957839 %d %d %d 1\n'
varname_ssh = ' water level [m]\n'

fout_ssh = 'zinit.dat' 

ff 	= open(fout_ssh,'w')

nvars	= 1

ff.write(header % (nkn, 1,nvars))
ff.write(date_init_str)
ff.write(varname_ssh)
for nn in range(nkn):
        ff.write('%.8f\n' % ssh[nn])
ff.close()

