import numpy as np
from pandas.tseries.offsets import DateOffset

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'

# df=xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/20000101.ocean_month_2000_09.nc',decode_times=False)


ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[0:12]
for itr in range(2,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2,decode_times=False)

# Fram Strait indices are j=724:844, i = 674
# monthly transport in Sv
Tr_Fr = df.umo[:,:,724:844,1424].sum(('zl','yh'))*1e-9

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Tr_Br.shape[0])
Tr_Fr['Time'] = time
dfs = Tr_Fr.to_dataset(name='Fram_volume_tr')
dfs.to_netcdf('Fram_volume_transport.nc')
