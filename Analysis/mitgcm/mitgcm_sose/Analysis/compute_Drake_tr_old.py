import xmitgcm

root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'
project_name = 'Exp01_0'
outname = project_name + '_Drake_Tr.nc'
fnames = root_folder+project_name+'/*UVELMASS*'
list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)


# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
# list = sorted(glob.glob(fnames))
# df = xr.open_mfdataset(list)

voltr = df.UVELMASS*df.dyG*df.drF
# location of the Drake Passage
Tr = voltr[:,:,525:761,3555]


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
