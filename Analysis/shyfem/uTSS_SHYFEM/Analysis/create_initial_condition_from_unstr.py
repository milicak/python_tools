from netCDF4 import Dataset
import sys
import numpy as np
from datetime import datetime, timedelta
from shyfem_utils import extend_last_value


date_format     = '%Y%m%d%H%M%S'
# filename = '/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/uTSS_lobc_chunk_1460.nos.nc'
# what = 'nos' 
filename = '/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/uTSS_lobc_chunk_1460.ous.nc'
what = 'ous'
date0_string = '20200101000000'
time_level = int(0)

date0   = datetime.strptime(date0_string,date_format)

nc      = Dataset(filename)
time_seconds            = nc.variables['time'][:]
print ('time seconds', time_seconds)
print ('type', type(time_seconds))
time_selected   = time_seconds[time_level]
print ('time selected', time_selected)
print ('type2', type(time_selected))
nkn     = nc.dimensions['node'].size
print ('nkn',nkn)

levels  = nc.variables['level'][:]
if what == 'nos':
    print ('Opening nos file')
    temp        = nc.variables['temperature'][time_level,:]
    salt        = nc.variables['salinity'][time_level,:]
    #temp[temp == 0.0] = -999.0    
    #salt[salt == 0.0] = -999.0    
    temp = extend_last_value(temp)
    salt = extend_last_value(salt)
elif what == 'ous':    
    print ('Opening ous file')
    ssh         = nc.variables['water_level'][time_level,:]
    uvel        = nc.variables['u_velocity'][time_level,:]
    vvel        = nc.variables['v_velocity'][time_level,:]
    #uvel[uvel == 0.0] = -999.0    
    #vvel[vvel == 0.0] = -999.0    
    uvel = extend_last_value(uvel)
    vvel = extend_last_value(vvel)
else:
    print('Unrecognized type of file')
    sys.exit(1)

nc.close() 

# date1   = date0 + timedelta(seconds=time_selected)
date1   = date0 

# set in marmara velocities below threshold to zero
df = xr.open_dataset(filename)
mask = xr.where(((df.longitude > 26.0) & 
                 (df.longitude < 30) & 
                 (df.latitude > 40.2) &
                 (df.latitude < 41.2)),0,1) 
mask2=np.tile(mask,(23,1))  
uvel[:,70::] = uvel[:,70::]*np.transpose(mask2)
vvel[:,70::] = vvel[:,70::]*np.transpose(mask2)

### write to file
header          = '0 2 957839 %d %d %d 1\n'
bottom_levels_list = ' '.join(str(ii) for ii in levels) + '\n'
temp_varname    = ' temperature [C]\n'
salt_varname    = ' salinity [psu]\n'
ssh_varname     = ' water level [m]\n' 
uvel_varname    = ' u-velocity [m/s]\n'
vvel_varname    = ' v-velocity [m/s]\n'

date_format_out = '%Y%m%d %H%M%S\n'

nz              = len(levels)

ssh_filename    = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/boundin_%s.dat' % (date1.strftime('%Y%m%d_%H%M%S'))
vel_filename    = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/vel3din_%s.dat' % (date1.strftime('%Y%m%d_%H%M%S'))
temp_filename   = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/tempin_%s.dat' % (date1.strftime('%Y%m%d_%H%M%S'))
salt_filename   = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/saltin_%s.dat' % (date1.strftime('%Y%m%d_%H%M%S'))

print ('writing files...')
#### write file
if what == 'nos':
    ### TEMPERATURE
    ff  = open(temp_filename,'w')
    ff.write(header % (nkn, nz, 1))
    ff.write(date1.strftime(date_format_out))
    ff.write(bottom_levels_list)
    ff.write(temp_varname)
    string0 = '%s -999.0 ' % nz 
    for nn in range(nkn):
        stringg = string0 + ' '.join(' %.12f ' % ii for ii in temp[nn,:].tolist()) + '\n'
        ff.write(stringg)
    ff.close()
    ### SALINITY
    ff  = open(salt_filename,'w')
    ff.write(header % (nkn, nz, 1))
    ff.write(date1.strftime(date_format_out))
    ff.write(bottom_levels_list)
    ff.write(salt_varname)
    string0 = '%s -999.0 ' % nz 
    for nn in range(nkn):
        stringg = string0 + ' '.join(' %.12f ' % ii for ii in salt[nn,:].tolist()) + '\n'
        ff.write(stringg)
    ff.close()
else:
    ### WATER LEVEL
    ff  = open(ssh_filename,'w')
    ff.write(header % (nkn, 1,1))
    ff.write(date1.strftime(date_format_out))
    ff.write(ssh_varname)
    for nn in range(nkn):
        ff.write('%.12f\n' % ssh[nn])
    ff.close()
    ### VELOCITIES
    ff  = open(vel_filename,'w')
    ff.write(header % (nkn, nz,2))
    ff.write(date1.strftime(date_format_out))
    ff.write(bottom_levels_list)
    ff.write(uvel_varname)
    string0 = '%s -999.0 ' % nz 
    for nn in range(nkn):
        stringg = string0 + ' '.join(' %.12f ' % ii for ii in uvel[nn,:].tolist()) + '\n'
        ff.write(stringg)
    ff.write(vvel_varname)
    for nn in range(nkn):
        stringg = string0 + ' '.join(' %.12f ' % ii for ii in vvel[nn,:].tolist()) + '\n'
        ff.write(stringg)
    ff.close()


