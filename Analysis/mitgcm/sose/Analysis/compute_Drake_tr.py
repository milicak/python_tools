import xmitgcm

root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
print(project_name)
outname = project_name + '_Drake_Tr.nc'
fnames = root_folder+project_name+'/*UVELMASS*'
grname = root_folder+project_name+'/grid.nc'
list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)

gr = xr.open_dataset(grname)

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
# list = sorted(glob.glob(fnames))
# df = xr.open_mfdataset(list)

df = df.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
voltr = df.UVELMASS*gr.dyG
voltr = voltr*gr.drF
# location of the Drake Passage
Tr = voltr[:,:,525:761,3555]

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
