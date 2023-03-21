import numpy as np
from matplotlib.path import Path

ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/*nos.nc'))

lonrec = np.array([27.525,28.775,28.775,27.525,27.525])
latrec = np.array([40.525,40.525,40.88,40.88,40.525])
vertices = np.vstack((lonrec,latrec))
mpath = Path( np.transpose(vertices) )

ds = xr.open_dataset(ls1[0])
lonlat = np.transpose(np.vstack((np.copy(ds.longitude),np.copy(ds.latitude))))
mask_flat = mpath.contains_points(lonlat)


dd = []
for count, ind in enumerate(ls1):
    print(count)
    ds = xr.open_dataset(ind)
    dnm = ds.temperature[:,:,0]*mask_flat
    dnm = dnm.where(dnm!=0)
    dnm2 = dnm.dropna('node')
    df = dnm2.to_dataset(name='shyfem_sst')
    if count == 0:
        dd = df
    else:
        dd = xr.merge([dd,df])


ls2 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT/*nos.nc'))
for count, ind in enumerate(ls2):
    print(count)
    ds = xr.open_dataset(ind)
    dnm = ds.temperature[:,:,0]*mask_flat
    dnm = dnm.where(dnm!=0)
    dnm2 = dnm.dropna('node')
    df = dnm2.to_dataset(name='shyfem_sst')
    df['level'] = 1.0
    dd = xr.merge([dd,df])


dd.to_netcdf('/work/opa/mi19918/Projects/nemo/Marmara/Marmara_shyfem_SST_2016_2021.nc')
