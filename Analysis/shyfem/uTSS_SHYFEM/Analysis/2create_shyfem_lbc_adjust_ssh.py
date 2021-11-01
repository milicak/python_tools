import numpy as np
import xarray as xr
import scipy
from datetime import datetime, timedelta
from my_sea_over_land import seaoverland_2d,seaoverland_1d 
from shyfem_utils import write_record_in_file1d, write_record_in_file2d, write_record_in_file2d_vel

date0 = datetime(2020,1,1,0)
dates = [date0 + timedelta(days=d) for d in range(366)]
# dates = [date0 + timedelta(days=d) for d in range(31)]
itr = 0

out_file_ssh1  = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/boundn_L1_adj.dat'
out_file_ssh2  = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/boundn_L2_adj.dat'
out_file_ssh3  = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/boundn_L3_adj.dat'
varname_ssh   = ' water level [m]\n'
fout_ssh1  = open(out_file_ssh1,'w')
fout_ssh2  = open(out_file_ssh2,'w')
fout_ssh3  = open(out_file_ssh3,'w')

df1 = xr.open_dataset('uTSS_boundary_valuesL1.nc')
df2 = xr.open_dataset('uTSS_boundary_valuesL2.nc')
df3 = xr.open_dataset('uTSS_boundary_valuesL3.nc')


nbound1 = df1.zos.shape[1] 
nbound2 = df2.zos.shape[1] 
nbound3 = df3.zos.shape[1] 

boundn1 = np.copy(df1.zos)
boundn2 = np.copy(df2.zos)
boundn3 = np.copy(df3.zos)

# calculate average over period of L1 and L2 (weighted average)
w1,w2 = float(nbound1)/float(nbound1+nbound2),float(nbound2)/float(nbound1+nbound2)

l1l2_ave = w1 * np.average(boundn1) + w2 * np.average(boundn2)
l1l2_min = w1 * np.min(boundn1) + w2 * np.min(boundn2)


# remove average from boundaries
boundn1a = boundn1 - l1l2_ave
boundn2a = boundn2 - l1l2_ave
boundn3a = boundn3 - l1l2_ave
    

for date in dates:
    print(date)
    transect_ssh1 = boundn1a[itr, :]
    transect_ssh2 = boundn2a[itr, :]
    transect_ssh3 = boundn3a[itr, :]
    itr += 1
    ## write recod in files
    write_record_in_file1d(fout_ssh1,date,varname_ssh,nbound1,1,1,np.copy(transect_ssh1))
    write_record_in_file1d(fout_ssh2,date,varname_ssh,nbound2,1,1,np.copy(transect_ssh2))
    write_record_in_file1d(fout_ssh3,date,varname_ssh,nbound3,1,1,np.copy(transect_ssh3))



fout_ssh1.close()
fout_ssh2.close()
fout_ssh3.close()
