import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
project_name = 'Exp05_0s'
outname = project_name + '_stream_function.nc'
fnames = root_folder+project_name+'/'

dtime = 100 # 100 or 80
dtime = 80 # 100 or 80

df = xmitgcm.open_mdsdataset(fnames,
                             prefix='UVELMASS',geometry='sphericalpolar',
                            ref_date='2006-12-31 12:0:0',delta_t=dtime)

df = df.sel(time=slice('2007-01-01','2013-12-31'))

mask = df.maskC[0,:,:]
mask = np.multiply(mask,1)
mask = mask.rename({"XC": "XG"})

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
# list = sorted(glob.glob(fnames))
# df = xr.open_mfdataset(list)

Ubar = df.UVELMASS*df.drF
Ubar = Ubar.sum('Z')
dnm = Ubar*df.dyG
dnm
aa = dnm.cumsum('YC')
aa = aa.mean('time')
dnm = aa*mask.data
ds = dnm.to_dataset(name='stream_function')
ds.to_netcdf(outname)




