import numpy as np
import scipy.io

time = 17.5

df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus_oldlonlat/NEMO_TS_IC_v2.nc')

lon_rho = np.flipud(np.transpose(np.copy(df.lon)))
lat_rho = np.flipud(np.transpose(np.copy(df.lat)))
nk = df.temp.shape[1]

var1 = np.flipud(np.transpose(np.copy(df.temp[0,0,:,:])))
nj, ni = var1.shape

var = np.zeros((nk,nj,ni))
varsalt = np.zeros((nk,nj,ni))
for k in np.arange(0,nk):
    var[k,:,:] = np.flipud(np.transpose(np.copy(df.temp[0,k,:,:])))
    varsalt[k,:,:] = np.flipud(np.transpose(np.copy(df.salt[0,k,:,:])))


var =np.reshape(var,(1,nk,nj,ni))
varsalt =np.reshape(varsalt,(1,nk,nj,ni))

# Create a mosaic file
fout = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/' + 'NEMO_TS_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')
# Dimensions
rg.createDimension('time', None)
rg.createDimension('z_l',nk)
rg.createDimension('nxp',ni)
rg.createDimension('nyp',nj)
# Variables
hnx = rg.createVariable('nxp', 'int32', ('nxp',))
hny = rg.createVariable('nyp', 'int32', ('nyp',))
hx = rg.createVariable('lon','float32',('nyp','nxp',))
hx.units = 'degrees'
hy = rg.createVariable('lat','float32',('nyp','nxp',))
hy.units = 'degrees'
hz = rg.createVariable('z_l','float32',('z_l',))
hz.units = 'meters'
tempvar  = rg.createVariable('temp','float32',('time','z_l','nyp','nxp',))
tempvar.units = 'celcius'
tempvar.missing_val = 1e20
tempvar._FillValue = 1e20
saltvar  = rg.createVariable('salt','float32',('time','z_l','nyp','nxp',))
saltvar.units = 'psu'
saltvar.missing_val = 1e20
saltvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1958-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
hz[:] = np.copy(df.z_l)
tempvar[:] = var
saltvar[:] = varsalt
hnx[:] = np.arange(0,ni)
hny[:] = np.arange(0,nj)
htime = time
rg.close()


