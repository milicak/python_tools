import numpy as np
from pandas.tseries.offsets import DateOffset

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'

# df=xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/20000101.ocean_month_2000_09.nc',decode_times=False)


ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[0:12]
for itr in range(2,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2,decode_times=False)

# Davis Strait indices are j=345:500, i = 1450
# monthly transport in TW
Ht_Dv = (df.T_adx_2d[:,345:500,1450].sum('yh') +
      df.T_diffx_2d[:,345:500,1450].sum('yh') ) * 1e-12

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Ht_Dv.shape[0])
Ht_Dv['Time'] = time
dfs = Ht_Dv.to_dataset(name='Davis_heat_tr')
dfs.to_netcdf('Davis_heat_transport.nc')
