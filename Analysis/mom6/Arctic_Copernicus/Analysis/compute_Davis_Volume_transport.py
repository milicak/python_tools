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

# ue = np.cos(angle_dx.data*deg_rad)*urho+np.sin(angle_dx.data*deg_rad)*vrho
# vn = -np.sin(angle_dx.data*deg_rad)*urho+np.cos(angle_dx.data*deg_rad)*vrho

ue = np.cos(angle_dx.data*deg_rad)*urho-np.sin(angle_dx.data*deg_rad)*vrho
vn = np.sin(angle_dx.data*deg_rad)*urho+np.cos(angle_dx.data*deg_rad)*vrho

# Davis Strait
# monthly transport in Sv
aa = np.arange(1450,1501)
bb = np.arange(350,401)
y = xr.DataArray(bb, dims='points')
x = xr.DataArray(aa, dims='points')
# UTr = df.umo.isel(xq=x,yh=y)
VTr = vn.isel(xh=x,yh=y)
Tr_Dv = VTr.sum(('zl','points'))*1e-9

Tr_Dv = df.umo[:,:,345:500,1450].sum(('zl','yh'))*1e-9

# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=Tr_Dv.shape[0])
Tr_Dv['Time'] = time
dfs = Tr_Dv.to_dataset(name='Davis_volume_tr')
dfs.to_netcdf('Davis_volume_transport.nc')
