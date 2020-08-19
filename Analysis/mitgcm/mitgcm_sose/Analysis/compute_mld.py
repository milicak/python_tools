import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
project_name = 'Exp01_0s'
outname = project_name + '_taux.nc'
fnames = root_folder+project_name+'/'

dtime = 100 # 100 or 80
# dtime = 80 # 100 or 80

df = xmitgcm.open_mdsdataset(fnames,
                             prefix='KPPhbl',geometry='sphericalpolar',
                            ref_date='2006-12-31 12:0:0',delta_t=dtime)
df = xmitgcm.open_mdsdataset(fnames,
                             prefix='SIarea',geometry='sphericalpolar',
                            ref_date='2006-12-31 12:0:0',delta_t=dtime)

df = df.sel(time=slice('2007-01-01','2013-12-31'))
ds_by_season = df['KPPhbl'].groupby('time.season').mean('time')
ds = ds_by_season.sel(season='JJA')
ds = ds['mld']
ds = ds.to_dataset(name='mld')
ds.to_netcdf(outname)




