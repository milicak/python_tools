import numpy as np
import xarray as xr



# grab seawifs data from Udel thredds, subset to our domain, resample to make monthly climatology. Note that resampling will require some decent computation time.
# we only want Chlorophlyll, so we subset that variable as well

wifs = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/SEAWIFSDAILY9KM.nc").sel(lat=slice(-29,81),lon=slice(-104,57))
wifs = xr.open_dataset('https://basin.ceoe.udel.edu/thredds/dodsC/SEAWIFSMISSION9KM.nc')
wifs = wifs['chlor_a']
wifs_mon = wifs.groupby('time.month').mean('time')
# wifs.to_netcdf("/okyanus/users/milicak/dataset/MOM6/NA12/seawifs-clim-1997-2010.nc", format='NETCDF3_64BIT')
wifs.to_netcdf("/okyanus/users/milicak/dataset/MOM6/NA12/seawifs-clim-1997-2010.nc")

