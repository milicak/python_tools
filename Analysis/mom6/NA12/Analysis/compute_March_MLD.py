from xgcm import Grid
from mpl_toolkits.basemap import Basemap

gr = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_geometry.nc')
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})

ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/*ocean_month*.nc'))
keyword = 'snap'
for i in ls1[:]:
    if i.find(keyword.lower()) != -1:
        ls1.remove(i)

keyword = '_z_'
for i in ls1[:]:
    if i.find(keyword.lower()) != -1:
        ls1.remove(i)

# skip 1996 and 1997
ls1 = ls1[24:]
# get only march values
df = xr.open_mfdataset(ls1[2::12])
ds = df.MLD_003.mean('time')
ds1 = ds.to_dataset(name='MLD_003')
ds1.to_netcdf('MOM6_mean_MLD_003.nc')
