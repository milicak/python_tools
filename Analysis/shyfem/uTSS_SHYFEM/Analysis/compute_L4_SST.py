import numpy as np


root_folder = '/data/inputs/metocean/historical/obs/satellite/SST/CNR/L4/day/'


ls1 = sorted(glob.glob('/data/inputs/metocean/historical/obs/satellite/SST/CNR/L4/day/*/*/*UHR_NRT-MED-v02.0-fv02.0*'))


df = xr.open_mfdataset(ls1)
ds = df.analysed_sst[:,1200:1300,5400:5800] - 273.15
ds.to_netcdf('/work/opa/mi19918/Projects/nemo/Marmara/Marmara_L4_SST_2008_2022.nc')


# ls1 = sorted(glob.glob('/data/inputs/metocean/historical/obs/satellite/SST/CNR/L4/day/2016/*/*L4_GHRSST-SSTdepth-OSTIA-GLOB**'))
ls1 = sorted(glob.glob('/data/inputs/metocean/historical/obs/satellite/SST/CNR/L4/day/*/*/*GLOB**'))
df = xr.open_mfdataset(ls1)

dd = []
for count, ind in enumerate(ls1):
    print(count)
    df = xr.open_dataset(ind)
    ds = df.analysed_sst[:,2600:2623,4130:4200] - 273.15
    ds = ds.to_dataset(name='analysed_sst')
    if count == 0:
        dd = ds
    else:
        dd = xr.merge([dd,ds])



dd.to_netcdf('/work/opa/mi19918/Projects/nemo/Marmara/Marmara_global_L4_SST_2004_2020.nc')
