import numpy as np


root_folder = '/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/reanalysis/day/'

# gr = xr.open_dataset('/data/opa/bs-mod/experiments/bs-simu_6.6_zeus/36_1d_20140903_20140909_grid_T.nc')
gr = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/mesh_mask_Aug19_31-5-25.nc')
gr = gr.rename({'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'})
# year = 1993
# files = root_folder + np.str(year) + '/*/*RFVL*BSe3r1*'      
files = root_folder + '*/*/*TEMP*BSe3r1*'      
ls1 = sorted(glob.glob(files)) 
files = root_folder + '*/*/*PSAL*BSe3r1*'      
ls2 = sorted(glob.glob(files)) 
comp = dict(zlib=True, complevel=5)
# ind = 0


# 1993 - 1996
lst = ls1[0:365*4+1]  
lss = ls2[0:365*4+1]  
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  
skip = 10
plt.scatter(dfs.so[:,:,::10,::10],dft.thetao[:,:,::10,::10]);plt.xlim(15,25);plt.ylim(6,12);

# 1997 - 2000
lst = ls1[365*4+1:365*8+2]    
lss = ls2[365*4+1:365*8+2]    
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  

# 2001 - 2004
lst = ls1[365*8+2:365*12+3]    
lss = ls2[365*8+2:365*12+3]    
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  

# 2005 - 2008
lst = ls1[365*12+3:365*16+4] 
lss = ls2[365*12+3:365*16+4] 
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  


# 2009 - 2012
lst = ls1[365*16+4:365*20+5] 
lss = ls2[365*16+4:365*20+5] 
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  

# 2013 - 2016
lst = ls1[365*20+5:365*24+6] 
lss = ls2[365*20+5:365*24+6] 
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  

# 2017 - 2020
lst = ls1[365*24+6:365*28+6] 
lss = ls2[365*24+6:365*28+6] 
dft = xr.open_mfdataset(lst)
dfs = xr.open_mfdataset(lss)
dft = dft.isel(lat=95,lon=98) 
dfs = dfs.isel(lat=95,lon=98) 
dft['so']=dfs.so  


skip = 10
plt.scatter(dfs.so[:,:,::10,::10],dft.thetao[:,:,::10,::10]);plt.xlim(15,25);plt.ylim(6,12);
     



for ind in range(0,len(ls1)):
    dft = xr.open_dataset(ls1[ind])  
    dfs = xr.open_dataset(ls2[ind])  
    dft = dft.isel(lat=95,lon=98) 
    dfs = dfs.isel(lat=95,lon=98) 
    dft['so']=dfs.so  
    # dfs1 = dfs.isel(time=0) 
    # dft1 = dft.isel(time=0) 
    # # TS values transport                                           
    Trxreverse_zsum = Trxreverse_zsum.to_dataset(name='MOC_meridional')
    print(ls1[ind][79:87])
    outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/TS_diagram_43_5N31E_BSe3r1_' + ls1[ind][79:87] + '.nc'
    encoding = {var: comp for var in dft.data_vars}                             
    dft.to_netcdf(outname, encoding=encoding)                                     
