/okyanus/users/milicak/dataset/soda3.7.2shrt/ocn/ORIGINAL/ocean
list=sorted(glob.glob('*2007*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2007_temp_northernBC.nc')

list=sorted(glob.glob('*2008*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2008_temp_northernBC.nc')

list=sorted(glob.glob('*2009*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2009_temp_northernBC.nc')

list=sorted(glob.glob('*2010*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2010_temp_northernBC.nc')

list=sorted(glob.glob('*2011*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2011_temp_northernBC.nc')

list=sorted(glob.glob('*2012*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2012_temp_northernBC.nc')

list=sorted(glob.glob('*2013*'))
ds=xr.open_mfdataset(list)['temp']
ds2=ds[:,:,388,:]
del ds2.attrs['coordinates']
ds2.to_netcdf('2013_temp_northernBC.nc')

