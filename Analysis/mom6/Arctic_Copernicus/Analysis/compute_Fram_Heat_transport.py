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
# Fram Strait indices are j=724:844, i = 674
# monthly transport in TW
Ht_Fr = (df.T_adx_2d[:,724:844,1424].sum('yh') +
      df.T_diffx_2d[:,724:844,1424].sum('yh') ) * 1e-12

Ht_Fr_Aw = (df.T_adx_2d[:,724:844,1424].where(df.T_adx_2d[:,724:844,1424]>0).sum('yh')
            +df.T_diffx_2d[:,724:844,1424].where(df.T_diffx_2d[:,724:844,1424]>0).sum('yh')
            ) * 1e-12

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Ht_Fr.shape[0])
Ht_Fr['Time'] = time
dfs = Ht_Fr.to_dataset(name='Fram_heat_tr')
dfs.to_netcdf('Fram_heat_transport.nc')


Ht_Fr_Aw['Time'] = time
dfs = Ht_Fr_Aw.to_dataset(name='Fram_AW_heat_tr')
dfs.to_netcdf('Fram_AW_heat_transport.nc')
