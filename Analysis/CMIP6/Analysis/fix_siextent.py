import xarray as xr
import numpy as np
import glob
from datetime import datetime

#from dask.distributed import Client, LocalCluster
#
from OMIP_utils import annual_mean

if __name__ == '__main__':
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=8, memory_limit=25e9, local_dir='/cluster/projects/nn2345k/documentation_figures_BLOM/diagnostics_coupled/NorESM2-LM/output/')
    client = Client(cluster)


#
outpath = '/cluster/projects/nn2345k/documentation_figures_BLOM/diagnostics_coupled/NorESM2-LM/output/'
datestr = '20191209_fixed'

#cluster = LocalCluster(n_workers=6, memory_limit=15e9, local_dir=outpath)
#client  = Client(cluster)
#
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

print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
#
siextentn = (ice.tarea).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#
for var in [siextentn,siextents]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+datestr+'.nc')

ice.close()
print('sea ice done: ' + expID)


#-------------------------------------------
# for the NorESM2-LM historical1
expID='historical1'
case1='NHIST_f19_tn14_20190625/'
case2='NHIST_f19_tn14_20190710/'

fnames_ocn1 = sorted(glob.glob(fpath+case1+'ocn/hist/*.micom.hm*.nc'))
fnames_ice1 = sorted(glob.glob(fpath+case1+'ice/hist/*cice.h.????-??.nc'))
fnames_ocn2 = sorted(glob.glob(fpath+case2+'ocn/hist/*.micom.hm*.nc'))
fnames_ice2 = sorted(glob.glob(fpath+case2+'ice/hist/*cice.h.????-??.nc'))

fnames_ocn=fnames_ocn1+fnames_ocn2
fnames_ice=fnames_ice1+fnames_ice2

year_range1=slice('1850','2014') # year range for the whole historical period
year_range2=slice('1980','2009') # year range for the certain historical period


print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
siextentn = (ice.tarea).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#
#
for var in [siextentn,siextents]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+datestr+'.nc')

ice.close()
print('sea ice done: ' +expID)


#-------------------------------------------
# for the NorESM2-LM historical2
expID='historical2'
case1='NHIST_02_f19_tn14_20190801/'
case2='NHIST_02_f19_tn14_20190813/'

fnames_ocn1 = sorted(glob.glob(fpath+case1+'ocn/hist/*.micom.hm*.nc'))
fnames_ice1 = sorted(glob.glob(fpath+case1+'ice/hist/*cice.h.????-??.nc'))
fnames_ocn2 = sorted(glob.glob(fpath+case2+'ocn/hist/*.micom.hm*.nc'))
fnames_ice2 = sorted(glob.glob(fpath+case2+'ice/hist/*cice.h.????-??.nc'))

fnames_ocn=fnames_ocn1+fnames_ocn2
fnames_ice=fnames_ice1+fnames_ice2

year_range1=slice('1850','2014') # year range for the whole historical period
year_range2=slice('1980','2009') # year range for the certain historical period


print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
#
siextentn = (ice.tarea).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#
for var in [siextentn,siextents]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+datestr+'.nc')

ice.close()
print('sea ice done: ' +expID)

#-------------------------------------------
# for the NorESM2-LM historical3
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


print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
siextentn = (ice.tarea).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#
#
for var in [siextentn,siextents]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+datestr+'.nc')

ice.close()
print('sea ice done: ' +expID)

#-------------------------------------------
# for the NorESM2-MM historical1
expID='historical1'
case1='NHISTfrc2_f09_tn14_20191001/'
case2='NHISTfrc2_f09_tn14_20191025/'

#fnames_ocn1 = sorted(glob.glob(fpath+case1+'ocn/hist/*.micom.hm*.nc'))
fnames_ice1 = sorted(glob.glob(fpath+case1+'ice/hist/*cice.h.????-??.nc'))
#fnames_ocn2 = sorted(glob.glob(fpath+case2+'ocn/hist/*.micom.hm*.nc'))
#fnames_ice2 = sorted(glob.glob(fpath+case2+'ice/hist/*cice.h.????-??.nc'))

#fnames_ocn=fnames_ocn1+fnames_ocn2
#fnames_ice=fnames_ice1+fnames_ice2
fnames_ice=fnames_ice1

year_range1=slice('1850','2014') # year range for the whole historical period
year_range2=slice('1980','2009') # year range for the certain historical period


print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
siextentn = (ice.tarea).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#
#
for var in [siextentn,siextents]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+'NHISTfrc2_f09_tn14_20191001_'+datestr+'.nc')

ice.close()
print('sea ice done: ' +expID)


#-------------------------------------------
# for the NorESM2-MM historical1 - part II
expID='historical1'
case1='NHISTfrc2_f09_tn14_20191001/'
case2='NHISTfrc2_f09_tn14_20191025/'

#fnames_ocn1 = sorted(glob.glob(fpath+case1+'ocn/hist/*.micom.hm*.nc'))
#fnames_ice1 = sorted(glob.glob(fpath+case1+'ice/hist/*cice.h.????-??.nc'))
#fnames_ocn2 = sorted(glob.glob(fpath+case2+'ocn/hist/*.micom.hm*.nc'))
fpath='/cluster/NS9560K/noresm/cases/'
fnames_ice2 = sorted(glob.glob(fpath+case2+'ice/hist/*cice.h.????-??.nc'))

#fnames_ocn=fnames_ocn1+fnames_ocn2
#fnames_ice=fnames_ice1+fnames_ice2
fnames_ice=fnames_ice2

year_range1=slice('1850','2014') # year range for the whole historical period
year_range2=slice('1980','2009') # year range for the certain historical period


print('load sea ice')
ice = xr.open_mfdataset(fnames_ice, combine='nested', concat_dim='time', parallel=True)

#cguo
newtime = ice.time_bounds.isel(d2=1).values.astype(datetime)-0.5*(ice.time_bounds.isel(d2=1).values.astype(datetime)-ice.time_bounds.isel(d2=0).values.astype(datetime))
ice=ice.assign_coords(time=newtime)

#
siextentn = (ice.tarea).sel(time=year_range1).where(ice.TLAT>0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextentn')
siextents = (ice.tarea).sel(time=year_range1).where(ice.TLAT<0).where(ice.aice.sel(time=year_range1)>.15).sum(dim=('ni','nj')).rename('siextents')
#
#
#
for var in [siextentn,siextents]:
    var.to_dataset(name=var.name).to_netcdf(outpath+var.name+'_'+expID+'_'+'NHISTfrc2_f09_tn14_20191025_'+datestr+'.nc')

ice.close()
print('sea ice done: ' +expID)


