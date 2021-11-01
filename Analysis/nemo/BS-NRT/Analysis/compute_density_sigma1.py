import numpy as np
import gsw


# from 1993 - 2019
root_folder = '/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/reanalysis/day/'
files = root_folder + '*/*/*TEMP*BSe3r1*'      
ls1 = sorted(glob.glob(files)) 
files = root_folder + '*/*/*PSAL*BSe3r1*'      
ls2 = sorted(glob.glob(files)) 
# for 2020
root_folder = '/data/products/BSFS/bs-rea_v3.1/cmems/reanalysis_daily_mean/'
files = root_folder + '2020/*TEMP*BSe3r1*'      
ls1 = sorted(glob.glob(files)) 
files = root_folder + '2020/*PSAL*BSe3r1*'      
ls2 = sorted(glob.glob(files)) 

# gr = xr.open_dataset('/data/opa/bs-mod/experiments/bs-simu_6.6_zeus/36_1d_20140903_20140909_grid_T.nc')
gr = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/mesh_mask_Aug19_31-5-25.nc')
gr = gr.rename({'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'})

comp = dict(zlib=True, complevel=5)

for ind in range(0,len(ls1)):
    dft = xr.open_dataset(ls1[ind])  
    dfs = xr.open_dataset(ls2[ind])  
    dfs1 = dfs.isel(time=0) 
    dft1 = dft.isel(time=0) 
    sigma1 = xr.apply_ufunc(gsw.sigma1, dfs1.so, dft1.thetao,                  
                            dask='parallelized', output_dtypes=[dfs1.so.dtype])
    df1 = sigma1.to_dataset(name='sigma1')                                      
    # from 1993 - 2019
    print(ls1[ind][79:87])
    outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/tmp/sigma1_BSe3r1_' + ls1[ind][79:87] + '.nc'
    # for 2020
    # print(ls1[ind][65:73])
    # outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/tmp/sigma1_BSe3r1_' + ls1[ind][65:73] + '.nc'
    encoding = {var: comp for var in df1.data_vars}                             
    df1.to_netcdf(outname, encoding=encoding)                                     

