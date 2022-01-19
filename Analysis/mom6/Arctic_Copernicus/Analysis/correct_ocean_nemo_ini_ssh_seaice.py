import numpy as np
import scipy.io

time = 17.5

df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus_oldlonlat/NEMO_ssh_IC.nc')
ds = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus_oldlonlat/NEMO_seaice_IC.nc')

lon_rho = np.flipud(np.transpose(np.copy(df.lon)))
lat_rho = np.flipud(np.transpose(np.copy(df.lat)))
ssh = np.flipud(np.transpose(np.copy(df.ssh[0,:,:])))
nj, ni = ssh.shape
hice = np.flipud(np.transpose(np.copy(ds.hice[0,:,:])))
aice = np.flipud(np.transpose(np.copy(ds.aice[0,:,:])))

ssh =np.reshape(ssh,(1,nj,ni))
aice =np.reshape(aice,(1,nj,ni))
hice =np.reshape(hice,(1,nj,ni))

# Create a mosaic file
fout = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/' + 'NEMO_ssh_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')
# Dimensions
rg.createDimension('time', None)
rg.createDimension('nxp',ni)
rg.createDimension('nyp',nj)
# Variables
hx = rg.createVariable('lon','float32',('nyp','nxp',))
hx.units = 'degrees'
hy = rg.createVariable('lat','float32',('nyp','nxp',))
hy.units = 'degrees'
sshvar  = rg.createVariable('ssh','float32',('time','nyp','nxp',))
sshvar.units = 'meters'
sshvar.missing_val = 1e20
sshvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1958-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
sshvar[:] = ssh
htime = time
rg.close()


# Create a mosaic file
fout = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/' + 'NEMO_seaice_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')
# Dimensions
rg.createDimension('time', None)
rg.createDimension('nxp',ni)
rg.createDimension('nyp',nj)
# Variables
hx = rg.createVariable('lon','float32',('nyp','nxp',))
hx.units = 'degrees'
hy = rg.createVariable('lat','float32',('nyp','nxp',))
hy.units = 'degrees'
aicevar  = rg.createVariable('aice','float32',('time','nyp','nxp',))
aicevar.units = 'concentration 0-1'
aicevar.missing_val = 1e20
aicevar._FillValue = 1e20
hicevar  = rg.createVariable('hice','float32',('time','nyp','nxp',))
hicevar.units = 'meters'
hicevar.missing_val = 1e20
hicevar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1958-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
aicevar[:] = aice
hicevar[:] = hice
htime = time
rg.close()

