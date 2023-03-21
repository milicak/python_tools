import numpy as np

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'

ls1 = sorted(glob.glob(root_folder + '*ice_month*'))
gr = xr.open_dataset(root_folder + 'sea_ice_geometry.nc')
gr = gr.rename({'lath':'yT','lonh':'xT'})

df = xr.open_mfdataset(ls1, decode_times=False)
time = pd.date_range("1996-01-01", freq="M", periods=264)
df['Time'] = time

aice = df.aice.fillna(0)
# sea ice extent criteria
SI = xr.where(aice<0.15,0,1)
SIice = SI*gr.Ah
AIice = aice*gr.Ah
AIicesum = AIice.sum(('yT','xT'))
SIicesum = SIice.sum(('yT','xT'))


ds = SIicesum.to_dataset(name='SI_extent')
ds['SI_area'] = AIicesum

ds.to_netcdf('seaice_extent_area.nc')

