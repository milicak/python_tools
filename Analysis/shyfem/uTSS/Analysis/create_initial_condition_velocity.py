import numpy as np
from datetime import datetime
from datetime import timedelta  
from fem_tools import load_shyfem_variable, load_shyfem_coords
import sys

startTime = datetime.now()

### to be set
ousnc_path	= 'uTSS_lobc_chunk_0640.ous.nc'
varindex = 0

time = xr.open_dataset(ousnc_path)['time']

# date_init 	= datetime(2016,1,1,0,0,0)	
date_init 	= datetime(2017,10,3,0,0,0)	

date_init_str 	= date_init.strftime('%Y%m%d %H%M%S\n')

lons,lats,elems,levels,time = load_shyfem_coords(ousnc_path)

n_levels	= len(levels)
nkn	= len(lons)

uvel	= load_shyfem_variable(ousnc_path,'u_velocity')
vvel	= load_shyfem_variable(ousnc_path,'v_velocity')

# select first time value
uvel	= uvel[0,:,:]
vvel	= vvel[0,:,:]	

print uvel.shape
print vvel.shape

###########################
### WRITE TO FILE
###########################
print 'write velocity file'

header = '0 2 957839 %d %d %d 1\n'
bottom_levels_list = ' '.join(str(ii) for ii in levels) + '\n'
varname_uvel = ' u-velocity [m/s]\n'
varname_vvel = ' v-velocity [m/s]\n'

fout_u = 'uvelin.dat' 
fout_v = 'vvelin.dat' 

ff 	= open(fout_u,'w')
ff2 	= open(fout_v,'w')

nvars	= 1

ff.write(header % (nkn, n_levels,nvars))
ff.write(date_init_str)
ff.write(bottom_levels_list)
ff.write(varname_uvel)
ff2.write(header % (nkn, n_levels,nvars))
ff2.write(date_init_str)
ff2.write(bottom_levels_list)
ff2.write(varname_vvel)
string0 = '%s -999.0 ' % n_levels
for nn in range(nkn):
        stringg_u = string0 + ' '.join(' %.8f ' % ii for ii in uvel[nn,:].tolist()) + '\n'
        stringg_v = string0 + ' '.join(' %.8f ' % ii for ii in vvel[nn,:].tolist()) + '\n'
        ff.write(stringg_u)
        ff2.write(stringg_v)
ff.close()
ff2.close()



print 'total execution time: ', datetime.now() - startTime




