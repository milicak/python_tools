import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
project_name = 'Exp01_0s'
outname = project_name + '_MOC_depth.nc'
fnames = root_folder+project_name+'/'

dtime = 100 # 100 or 80
# dtime = 80 # exp4.0s and exp5.0s

df = xmitgcm.open_mdsdataset(fnames,
                             prefix='VVELMASS',geometry='sphericalpolar',
                            ref_date='2006-12-31 12:0:0',delta_t=dtime)

df = df.sel(time=slice('2007-01-01','2013-12-31'))

# list = sorted(glob.glob(fnames))
# df = xr.open_mfdataset(list)

voltrV = df.VVELMASS*df.dxG*df.drF

Trx = voltrV.sum(('XC'))
Trxmean = Trx.mean('time')
# reverse z coordinate for different cumsum
Trxmeanreverse_z = Trxmean.reindex(Z=Trxmean.Z[::-1])

Trxmeanreverse_zsum = Trxmeanreverse_z.cumsum('Z')
Trxmeanreverse_zsum = Trxmeanreverse_zsum.reindex(Z=Trxmeanreverse_zsum.Z[::-1])
Trxmeanreverse_zsum = -Trxmeanreverse_zsum


ds = Trxmeanreverse_zsum.to_dataset(name='MOC_depth')
ds.to_netcdf(outname)


