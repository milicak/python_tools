import numpy as np
from pandas.tseries.offsets import DateOffset
from xgcm import Grid

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'
deg_rad = 1.0

ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[0:12]
for itr in range(2,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2,decode_times=False)
gr = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_geometry.nc')
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})

gr2 = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_hgrid.nc')
angle_dx = gr2.angle_dx[1::2,1::2]

# merge data
df = df.merge(gr)

# 2D grid
grid = Grid(df, coords={'X': {'center': 'xh', 'outer': 'xq'},
                        'Y': {'center': 'yh', 'outer': 'yq'},
                         }, periodic=[])

urho = grid.interp(df.umo, 'X', boundary='fill')
vrho = grid.interp(df.vmo, 'Y', boundary='fill')

ue = np.cos(angle_dx.data*deg_rad)*urho-np.sin(angle_dx.data*deg_rad)*vrho
vn = np.sin(angle_dx.data*deg_rad)*urho+np.cos(angle_dx.data*deg_rad)*vrho

# Fram Strait
# monthly transport in Sv
aa = np.arange(880,950)
bb = np.arange(1480,1620)
y = xr.DataArray(aa, dims='points')
x = xr.DataArray(bb[::2], dims='points')
x2 = xr.DataArray(bb[1::2], dims='points')

VTr1 = vn.isel(xh=x,yh=y)
Tr_BaV1 = VTr1.sum(('zl','points'))*1e-9
VTr2 = vn.isel(xh=x2,yh=y)
Tr_BaV2 = VTr2.sum(('zl','points'))*1e-9

UTr1 = ue.isel(xh=x,yh=y)
Tr_BaU1 = UTr1.sum(('zl','points'))*1e-9
UTr2 = ue.isel(xh=x2,yh=y)
Tr_BaU2 = UTr2.sum(('zl','points'))*1e-9

Tr_Ba = df.vmo[:,:,910,1457:1687].sum(('zl','xh'))*1e-9

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Tr_BaU1.shape[0])

dfs = Tr_BaV.to_dataset(name='Barents_V_volume_tr')
dfs['Time'] = time
dfs.to_netcdf('Barents_volume_transport.nc')


dfs = Tr_BaV1.to_dataset(name='Barents_V_volume_tr1')
dfs['Barents_V_volume_tr2'] = Tr_BaV2
dfs['Barents_U_volume_tr1'] = Tr_BaU1
dfs['Barents_U_volume_tr2'] = Tr_BaU2
dfs['Time'] = time
dfs.to_netcdf('Barents_volume_transport.nc')

