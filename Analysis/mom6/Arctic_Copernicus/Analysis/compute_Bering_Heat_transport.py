import numpy as np
from pandas.tseries.offsets import DateOffset

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'

# df=xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/20000101.ocean_month_2000_09.nc',decode_times=False)


ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[0:12]
for itr in range(2,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2,decode_times=False)

# Bering Strait indices are j=689:713, i = 674
# monthly transport in TW
Ht_Br = (df.T_adx_2d[:,689:713,674].sum('yh') +
      df.T_diffx_2d[:,689:713,674].sum('yh') ) * 1e-12

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Ht_Br.shape[0])
Ht_Br['Time'] = time
dfs = Ht_Br.to_dataset(name='Bering_heat_tr')
dfs.to_netcdf('Bering_heat_transport.nc')
