import numpy as np
from pandas.tseries.offsets import DateOffset

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'

# df=xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/20000101.ocean_month_2000_09.nc',decode_times=False)


ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[0:12]
for itr in range(2,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2,decode_times=False)

# Barents Sea indices are j=910, i = 1457:1687
# monthly transport in TW
Ht_Ba = (df.T_ady_2d[:,910,1457:1687].sum('xh') +
      df.T_diffy_2d[:,910,1457:1687].sum('xh') ) * 1e-12

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Ht_Ba.shape[0])
Ht_Ba['Time'] = time
dfs = Ht_Ba.to_dataset(name='Barents_heat_tr')
dfs.to_netcdf('Barents_heat_transport.nc')
