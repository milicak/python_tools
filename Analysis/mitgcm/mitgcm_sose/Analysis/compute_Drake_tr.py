import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
project_name = 'Exp05_0s'
outname = project_name + '_Drake_Tr.nc'
fnames = root_folder+project_name+'/'

dtime = 100 # 100 or 80
dtime = 80 # 100 or 80

df = xmitgcm.open_mdsdataset(fnames,
                             prefix='UVELMASS',geometry='sphericalpolar',
                            ref_date='2006-12-31 12:0:0',delta_t=dtime)

df = df.sel(time=slice('2007-01-01','2013-12-31'))

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
# list = sorted(glob.glob(fnames))
# df = xr.open_mfdataset(list)

voltr = df.UVELMASS*df.dyG*df.drF
# location of the Drake Passage
Tr = voltr[:,:,160:278,3550]


# Tr = voltr.isel(XG=3527,YC=125:280)
# old section
# Tr = voltr[:,:,125:280,3527]
# new section
# Tr = voltr[:,:,160:278,3550]


Trtime = Tr.sum(('Z','YC'))
ds = Trtime.to_dataset(name='DrakeTr')
ds.to_netcdf(outname)



# Tr2 = voltr[:,:,95:330,3480]
# Tr2time = Tr2.sum(('Z','YC'))

# Ubar=df3.UVELMASS*df3.drF
# Ubar=Ubar.sum('Z')
# dnm=Ubar*df3.dyG
# dnm
# aa=dnm.cumsum('YC')
#
