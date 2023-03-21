import numpy as np
import numpy.ma as ma
import glob
import sys
import pandas as pd
import glob
import xarray as xr

# /data/inputs/metocean/historical/model/atmos/ECMWF/IFS_010/analysis/6h/netcdf/2021/11/20211128-ECMWF---AM0100-MEDATL-b20211129_an-fv11.00.nc

# root_folder = '/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_010/analysis/6h/netcdf/'
root_folder = '/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/'
ls1 = []
for year in range(2017,2023):
    ls1.extend(sorted(glob.glob(root_folder + str(year) + '/*/' +
                                '*ECMWF---AM0125-MEDATL-*fv10*.nc')))
    # ls1.extend(sorted(glob.glob(root_folder + str(year) + '/*/' + '*ECMWF---AM0100-MEDATL-*.nc')))

# df = xr.open_mfdataset(ls1)

# df = xr.open_dataset('/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_010/analysis/6h/netcdf/2021/11/20211128-ECMWF---AM0100-MEDATL-b20211129_an-fv11.00.nc')

lon1 = 29
lat1 = 41
lon2 = 29.1
lat2 = 41.2
lat3 = 41.3

# this is foe IFS_010
# iind1 = 1040 # for lon1
# jind1 = 310 # for lat1
# iind2 = 1041 # for lon2
# jind2 = 312 # for lat2
# iind3 = 1041 # for lon3
# jind3 = 313 # for lat3
# iind4 = 1040 # for lon4
# jind4 = 309 # for lat4
# this is foe IFS_0125
iind1 = 832 # for lon1
jind1 = 248 # for lat1
iind2 = 833 # for lon2
jind2 = 250 # for lat2
iind3 = 833 # for lon3
jind3 = 251 # for lat3
iind4 = 832 # for lon4
jind4 = 247 # for lat4

for ind,fname in enumerate(ls1):
    print(ind)
    df = xr.open_dataset(fname)
    # latitude needs to be reindexed 
    df2 = df.V10M.reindex(lat=list(reversed(df.lat)))
    V1 = df2.isel(lon=iind1,lat=jind1)
    V1d = V1.resample(time='1D').mean('time')
    V2 = df2.isel(lon=iind2,lat=jind2)
    V2d = V2.resample(time='1D').mean('time')
    V3 = df2.isel(lon=iind3,lat=jind3)
    V3d = V3.resample(time='1D').mean('time')
    V4 = df2.isel(lon=iind4,lat=jind4)
    V4d = V4.resample(time='1D').mean('time')
    # latitude needs to be reindexed 
    df2 = df.MSL.reindex(lat=list(reversed(df.lat)))
    P1 = df2.isel(lon=iind1,lat=jind1)
    P1d = P1.resample(time='1D').mean('time')
    P2 = df2.isel(lon=iind2,lat=jind2)
    P2d = P2.resample(time='1D').mean('time')
    P3 = df2.isel(lon=iind3,lat=jind3)
    P3d = P3.resample(time='1D').mean('time')
    P4 = df2.isel(lon=iind4,lat=jind4)
    P4d = P4.resample(time='1D').mean('time')
    if(ind==0):
        ds1 = xr.Dataset({'V1d':V1d, 'P1d':P1d})
        ds2 = xr.Dataset({'V2d':V2d, 'P2d':P2d})
        ds3 = xr.Dataset({'V3d':V3d, 'P3d':P3d})
        ds4 = xr.Dataset({'V4d':V4d, 'P4d':P4d})
    else:
        tmp = xr.Dataset({'V1d':V1d, 'P1d':P1d})
        ds1 = xr.concat((ds1,tmp),'time')
        tmp = xr.Dataset({'V2d':V2d, 'P2d':P2d})
        ds2 = xr.concat((ds2,tmp),'time')
        tmp = xr.Dataset({'V3d':V3d, 'P3d':P3d})
        ds3 = xr.concat((ds3,tmp),'time')
        tmp = xr.Dataset({'V4d':V4d, 'P4d':P4d})
        ds4 = xr.concat((ds4,tmp),'time')

ds1.to_netcdf('era_blocking_v1p1.nc')
ds2.to_netcdf('era_blocking_v2p2.nc')
ds3.to_netcdf('era_blocking_v3p3.nc')
ds4.to_netcdf('era_blocking_v4p4.nc')
for ind,fname in enumerate(ls1):
    print(ind)
    df = xr.open_dataset(fname)
    # latitude needs to be reindexed 
    df2 = df.V10M.reindex(lat=list(reversed(df.lat)))
    V1 = df2.isel(lon=slice(827,836),lat=slice(247,251))
    V1d = V1.resample(time='1D').mean('time')
    V1d = V1d.mean(('lat','lon'))
    # latitude needs to be reindexed 
    df2 = df.MSL.reindex(lat=list(reversed(df.lat)))
    P1 = df2.isel(lon=slice(827,836),lat=slice(247,251))
    P1d = P1.resample(time='1D').mean('time')
    P1d = P1d.mean(('lat','lon'))
    if(ind==0):
        ds1 = xr.Dataset({'V1d':V1d, 'P1d':P1d})
    else:
        tmp = xr.Dataset({'V1d':V1d, 'P1d':P1d})
        ds1 = xr.concat((ds1,tmp),'time')


ds1.to_netcdf('era_blocking_mean.nc')


root_folder2 = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/tmp2/'
ls3 = sorted(glob.glob(root_folder2+'*v_velocity*'))
root_folder4 = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/tmp3/'
ls4 = sorted(glob.glob(root_folder4+'*v_velocity*'))
ls4 = ls4[275:]
kind1 = 17
kind2 = 51
lensim = len(ls3+ls4)
vtop = np.zeros(lensim)
vbottom = np.zeros(lensim)
blocking = np.zeros(lensim)
for itr in range(0,lensim):
        print(itr)
        if(itr<len(ls4)):
            itrt = itr+275
            fname = root_folder4 + 'REG_' + str(itrt) + '_' + str(itrt+1) + '_0_v_velocity_0.0005.nc'
        else:
            itrt = itr-len(ls4)
            fname = root_folder2 + 'REG_' + str(itrt) + '_' + str(itrt+1) + '_0_v_velocity_0.0005.nc'

        df = xr.open_dataset(fname)
        vtop[itr] = np.copy(df.v_velocity[0,kind1,570,436])
        vbottom[itr] = np.copy(df.v_velocity[0,kind2,570,436])
        # a criteris
        tmptop = vtop[itr] 
        tmpbottom = vbottom[itr] 
        if tmptop>-0.1:tmptop=1;
        if tmpbottom<0.1:tmpbottom=-1;
        tmp = np.sign(vtop[itr]*vbottom[itr])
        # tmp = np.sign(tmptop*tmpbottom)
        if tmp == -1:
                blocking[itr] = 0
        else:
                if vtop[itr] <0:
                        blocking[itr] = -1
                else:
                        blocking[itr] = 1



# save the blocking to netcdf
time = pd.date_range("2017-01-01", freq="D", periods=2079)
ds = xr.Dataset({"blocking": ("time", blocking), "time": time})


# blocking = blocking[:lensim]

# plt.plot(P2d[0:lensim]*np.abs(blocking),V2d[0:lensim]*np.abs(blocking),'*')
# plt.plot(blocking*np.abs(V2d[:lensim]))
#
# totalblocs = np.nansum(np.abs(blocking))
# mslb = np.zeros(totalblocs)
# vvelb = np.zeros(totalblocs)
# blockb = np.zeros(totalblocs)
# kind = 0
# for itr in range(0,lensim):
#     print(itr)
#     if blocking[itr]==1:
#         mslb[kind] = 
#
