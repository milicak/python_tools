import numpy as np
from datetime import datetime
from datetime import timedelta  
from fem_tools import load_shyfem_variable, load_shyfem_coords
import sys

startTime = datetime.now()

### to be set
nosnc_path 	= 'uTSS_lobc_chunk_0640.nos.nc'
ousnc_path	= 'uTSS_lobc_chunk_0640.ous.nc'
varindex = 0

time = xr.open_dataset(nosnc_path)['time']

# date_init 	= datetime(2016,1,1,0,0,0)	
date_init 	= datetime(2017,10,3,0,0,0)	

date_init_str 	= date_init.strftime('%Y%m%d %H%M%S\n')

lons,lats,elems,levels,time = load_shyfem_coords(nosnc_path)

n_levels	= len(levels)
nkn	= len(lons)

temp	= load_shyfem_variable(nosnc_path,'temperature')
salt	= load_shyfem_variable(nosnc_path,'salinity')

ssh	= load_shyfem_variable(ousnc_path,'water_level')
uvel	= load_shyfem_variable(ousnc_path,'u_velocity')
vvel	= load_shyfem_variable(ousnc_path,'v_velocity')

# select first time value
temp	= temp[0,:,:]
salt	= salt[0,:,:]
ssh	= ssh[0,:]
uvel	= uvel[0,:,:]
vvel	= vvel[0,:,:]	

print temp.shape
print salt.shape
print ssh.shape
print uvel.shape
print vvel.shape

###########################
### WRITE TO FILE
###########################


print 'writing active tracers files...'
### write to file
header = '0 2 957839 %d %d %d 1\n'
bottom_levels_list = ' '.join(str(ii) for ii in levels) + '\n'
varname_temp = 'temperature [C]\n'
varname_salt = 'salinity [psu]\n'

fout_temp = 'tempin.dat' 
fout_salt = 'saltin.dat' 

ff 	= open(fout_temp,'w')
ff2 	= open(fout_salt,'w')

nvars	= 1

ff.write(header % (nkn, n_levels,nvars))
ff.write(date_init_str)
ff.write(bottom_levels_list)
ff.write(varname_temp)
ff2.write(header % (nkn, n_levels,nvars))
ff2.write(date_init_str)
ff2.write(bottom_levels_list)
ff2.write(varname_salt)
string0 = '%s -999.0 ' % n_levels
for nn in range(nkn):
        stringg_temp = string0 + ' '.join(' %.8f ' % ii for ii in temp[nn,:].tolist()) + '\n'
        stringg_salt = string0 + ' '.join(' %.8f ' % ii for ii in salt[nn,:].tolist()) + '\n'
        ff.write(stringg_temp)
        ff2.write(stringg_salt)
ff.close()
ff2.close()

########
print 'write water level file..'
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

#########
print 'write velocity file'

header = '0 2 957839 %d %d %d 1\n'
bottom_levels_list = ' '.join(str(ii) for ii in levels) + '\n'
varname_uvel = ' u-velocity [m/s]\n'
varname_vvel = ' v-velocity [m/s]\n'

fout_uv = 'uvin.dat' 

ff 	= open(fout_uv,'w')

nvars	= 1

ff.write(header % (nkn, n_levels,nvars))
ff.write(date_init_str)
ff.write(bottom_levels_list)
ff.write(varname_uvel)
string0 = '%s -999.0 ' % n_levels
for nn in range(nkn):
        stringg_u = string0 + ' '.join(' %.8f ' % ii for ii in uvel[nn,:].tolist()) + '\n'
        ff.write(stringg_u)
ff.write(varname_vvel)
for nn in range(nkn):
        stringg_v = string0 + ' '.join(' %.5f ' % ii for ii in vvel[nn,:].tolist()) + '\n'
        ff.write(stringg_v)
ff.close()



print 'total execution time: ', datetime.now() - startTime




