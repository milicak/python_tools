import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/Blacksea_lonlat/'
project_name = 'JRAexp_01'
fnames = root_folder+project_name+'/'

outnameu = project_name + '_MOC_zonal_depth.nc'
outnamev = project_name + '_MOC_meridional_depth.nc'

dtime = 100 # 100 or 80

dfu = xmitgcm.open_mdsdataset(fnames,
                             prefix='UVELMASS',geometry='curvilinear',
                            ref_date='2009-12-31 12:0:0',delta_t=dtime)

dfv = xmitgcm.open_mdsdataset(fnames,
                             prefix='VVELMASS',geometry='curvilinear',
                            ref_date='2009-12-31 12:0:0',delta_t=dtime)

dfu = dfu.sel(time=slice('2010-01-01','2016-12-31'))
dfv = dfv.sel(time=slice('2010-01-01','2016-12-31'))

# voltrV = df.VVELMASS*df.dxG*df.drF
# Trx = voltrV.sum(('XC'))
# Trxsum = Trx.cumsum('Z')
# Trxsummean = Trxsum.mean('time')

# zonal transport
voltrU = dfu.UVELMASS*dfu.dyG*dfu.drF
Try = voltrU.sum(('j'))
Trymean = Try.mean('time')
# reverse z coordinate for different cumsum
Trymeanreverse_z = Trymean.reindex(k=Trymean.k[::-1])

Trymeansum = Trymean.cumsum('k')
Trymeanreverse_zsum = Trymeanreverse_z.cumsum('k')
Trymeanreverse_zsum = Trymeanreverse_zsum.reindex(k=Trymeanreverse_zsum.k[::-1])

# meridional transport
voltrV = dfv.VVELMASS*dfv.dxG*dfv.drF
Trx = voltrV.sum(('i'))
Trxmean = Trx.mean('time')
# reverse z coordinate for different cumsum
Trxmeanreverse_z = Trxmean.reindex(k=Trxmean.k[::-1])

Trxmeansum = Trxmean.cumsum('k')
Trxmeanreverse_zsum = Trxmeanreverse_z.cumsum('k')
Trxmeanreverse_zsum = Trxmeanreverse_zsum.reindex(k=Trxmeanreverse_zsum.k[::-1])

Trxmeanreverse_zsum = -Trxmeanreverse_zsum
Trymeanreverse_zsum = -Trymeanreverse_zsum

ds = Trymeanreverse_zsum.to_dataset(name='MOC_zonal_depth')
ds.to_netcdf(outnameu)

ds = Trxmeanreverse_zsum.to_dataset(name='MOC_meridional_depth')
ds.to_netcdf(outnamev)


