import xarray as xr
import numpy as np
import glob
from datetime import datetime

from dask.distributed import Client, LocalCluster
#
from OMIP_utils import annual_mean, annual_mean_loop, thermosteric_sealvl, thermosteric_sealvl_ref_period, upper_700m_mean

if __name__ == '__main__':
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=6, memory_limit=25e9, local_dir='/cluster/projects/nn2345k/documentation_figures_BLOM/diagnostics_coupled/NorESM2-LM/output/')
    client = Client(cluster)


#

#
outpath = '/cluster/projects/nn2345k/documentation_figures_BLOM/diagnostics_coupled/NorESM2-LM/output/'
datestr = '20191119'

#cluster = LocalCluster(n_workers=6, memory_limit=15e9, local_dir=outpath)
#client  = Client(cluster)
#
# LOAD GRID AND MASK FILES
gridpath = '/cluster/shared/noresm/inputdata/ocn/micom/tnx1v4/20170601/'
grid = xr.open_dataset(gridpath+'grid.nc')
mask = xr.open_dataset(gridpath+'ocean_regions.nc')

#
# LOAD DATA (LAZY)
fpath = '/tos-project3/NS9560K/noresm/cases/'

#
#-------------------------------------------
# for the NorESM2-LM piControl
expID='piControl'
case1='N1850_f19_tn14_20190621/'
case2='N1850_f19_tn14_20190722/'
case3='N1850_f19_tn14_20190802/'

fnames_ocn1 = sorted(glob.glob(fpath+case1+'ocn/hist/*.micom.hm*.nc'))
fnames_ice1 = sorted(glob.glob(fpath+case1+'ice/hist/*cice.h.????-??.nc'))
fnames_ocn2 = sorted(glob.glob(fpath+case2+'ocn/hist/*.micom.hm*.nc'))
fnames_ice2 = sorted(glob.glob(fpath+case2+'ice/hist/*cice.h.????-??.nc'))
fnames_ocn3 = sorted(glob.glob(fpath+case3+'ocn/hist/*.micom.hm*.nc'))
fnames_ice3 = sorted(glob.glob(fpath+case3+'ice/hist/*cice.h.????-??.nc'))

fnames_ocn=fnames_ocn1+fnames_ocn2+fnames_ocn3
fnames_ice=fnames_ice1+fnames_ice2+fnames_ice3

year_range1=slice('1601','2100') # year range for the whole 500 years period of piControl
year_range2=slice('2001','2100') # year range for the last 100 years of piControl


#-------------------------------------------
# for the NorESM2-LM historical
expID='historical1'
case1='NHIST_f19_tn14_20190625/'
case2='NHIST_f19_tn14_20190710/'

expID='historical2'
case1='NHIST_02_f19_tn14_20190801/'
case2='NHIST_02_f19_tn14_20190813/'

expID='historical3'
case1='NHIST_03_f19_tn14_20190801/'
case2='NHIST_03_f19_tn14_20190813/'

fnames_ocn1 = sorted(glob.glob(fpath+case1+'ocn/hist/*.micom.hm*.nc'))
fnames_ice1 = sorted(glob.glob(fpath+case1+'ice/hist/*cice.h.????-??.nc'))
fnames_ocn2 = sorted(glob.glob(fpath+case2+'ocn/hist/*.micom.hm*.nc'))
fnames_ice2 = sorted(glob.glob(fpath+case2+'ice/hist/*cice.h.????-??.nc'))

fnames_ocn=fnames_ocn1+fnames_ocn2
fnames_ice=fnames_ice1+fnames_ice2

year_range1=slice('1850','2014') # year range for the whole historical period
year_range2=slice('1980','2009') # year range for the certain historical period



#
print('load ocean data')
ocn = xr.open_mfdataset(fnames_ocn, combine='nested', concat_dim='time', parallel=True)
#
# ##################################
#          OCEAN FIELDS
# ##################################
#
# The following are fast to calculate
# just renaming and writing out data
#
# 1d timeseries
#
print('timeseries and surface fields')
thetaoga = annual_mean(ocn.tempga.sel(time=year_range1)).rename('thetaoga')
soga     = annual_mean(ocn.salnga.sel(time=year_range1)).rename('soga')
tosga    = annual_mean(ocn.sstga.sel(time=year_range1)).rename('tosga')
sosga    = annual_mean(ocn.sssga.sel(time=year_range1)).rename('sosga')
#
mfo      = annual_mean(ocn.voltr.sel(time=year_range1)).rename('mfo') #this includes all the passages 

# 2D timeseries
# z-y
amoc_rapid = annual_mean(ocn.mmflxd.isel(region=0).sel(lat=26).sel(time=year_range1).max(dim='depth'))
amoc_rapid = amoc_rapid.rename('amoc_rapid')
#
# x-y
# cguo
tos    = ocn.sst.sel(time=year_range2).mean(dim='time').rename('tos')
sos    = ocn.sss.sel(time=year_range2).mean(dim='time').rename('sos')


#tos_gn    = ocn.sst.isel(time=last_cycle).rename('tos_gn')
#tos       = tos_gn.isel(y=slice(0,-1)).rename('tos')
#sos_gn    = ocn.sss.isel(time=last_cycle).rename('sos_gn')
#sos       = sos_gn.isel(y=slice(0,-1)).rename('sos')
#zos_gn    = ocn.sealv.isel(time=last_cycle).rename('zos_gn')
#zos       = zos_gn.isel(y=slice(0,-1)).rename('zos')
#mlotst_gn = ocn.mlts.isel(time=last_cycle).rename('mlotst_gn')
#mlotst    = mlotst_gn.isel(y=slice(0,-1)).rename('mlotst')
#
# ###################
#   WRITE TO NETCDF
#for var in [thetaoga,soga,tosga,sosga,mfo,amoc_rapid]:
#for var in [amoc_rapid,tos,sos]:
for var in [mfo]:
    print(var)
    var.to_dataset().to_netcdf(outpath+var.name+'_'+expID+'_'+datestr+'.nc')

del thetaoga,soga,tosga,sosga,mfo,hfbasin,amoc_rapid,tos,sos,zos,mlotst
#


 
# GLOBAL MEAN T-S PROFILES
# define
mask3D      = ocn.templvl.isel(time=0).drop('time').notnull() # create a 3D mask
depth_bnds1 = xr.ones_like(ocn.templvl.isel(time=0).drop('time'))*ocn.depth_bnds.isel(bounds=1).isel(time=0).drop('time')
depth_bnds0 = xr.ones_like(ocn.templvl.isel(time=0).drop('time'))*ocn.depth_bnds.isel(bounds=0).isel(time=0).drop('time')
pdepth3D    = xr.ones_like(ocn.templvl.isel(time=0).drop('time'))*grid.pdepth
# 
dplvl       = (xr.concat([depth_bnds1,pdepth3D],dim='dum').min(dim='dum') - xr.concat([depth_bnds0,pdepth3D],dim='dum').min(dim='dum'))*ocn.pbot/grid.pdepth
thetao      = (ocn.templvl*mask3D*dplvl*grid.parea).sum(dim=('x','y'))/(dplvl*mask3D*grid.parea).sum(dim=('x','y'))
so          = (ocn.salnlvl*mask3D*dplvl*grid.parea).sum(dim=('x','y'))/(dplvl*mask3D*grid.parea).sum(dim=('x','y'))
# write out
thetao      = annual_mean_loop(thetao.rename('thetao'), dum_folder='/cluster/work/users/cgu025/tmp/thetao/', savepath=outpath+'thetao_'+datestr+'.nc')
del thetao
so          = annual_mean_loop(so.rename('so'), dum_folder='/cluster/work/users/cgu025/tmp/so/', savepath=outpath+'so_'+datestr+'.nc')
del so
#
ocn.close()
del ocn
print('ocn done')
#
# ################################
#         SEA ICE FIELDS
# ################################ 
#
print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
print('compute sea ice')
siarean = (ice.tarea*ice.aice).sel(time=year_range1).where(ice.TLAT>0).sum(dim=('ni','nj')).rename('siarean')
siareas = (ice.tarea*ice.aice).sel(time=year_range1).where(ice.TLAT<0).sum(dim=('ni','nj')).rename('siareas')
#
siextentn = (ice.tarea*ice.aice).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea*ice.aice).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#sivoln = annual_mean((ice.tarea*ice.hi).where(ice.TLAT>0).sum(dim=('ni','nj'))).rename('sivoln')
#sivols = annual_mean((ice.tarea*ice.hi).where(ice.TLAT<0).sum(dim=('ni','nj'))).rename('sivols')
#
siconc = ice.aice.sel(time=year_range2).groupby('time.month').mean(dim='time').rename('siconc')
sithick = ice.sithick.sel(time=year_range2).groupby('time.month').mean(dim='time').rename('sithick')

#
for var in [siarean,siareas,siextentn,siextents,siconc,sithick]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+datestr+'.nc')

ice.close()
print('sea ice done')
