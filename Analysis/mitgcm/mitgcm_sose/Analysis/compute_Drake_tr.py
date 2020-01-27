


root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'

project_name = 'Exp01_0'

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
fnames = root_folder+project_name+'/'+'*UVELMASS*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)

# location of the Drake Passage
iind = 3527

voltr = df.UVELMASS*df.dyG*df.drF
# Tr = voltr.isel(XG=3527,YC=125:280)
Tr = voltr[:,:,125:280,3527]
Trtime = Tr.sum(('Z','YC'))


Tr2 = voltr[:,:,95:330,3480]
Tr2time = Tr2.sum(('Z','YC'))

Ubar=df3.UVELMASS*df3.drF
Ubar=Ubar.sum('Z')
dnm=Ubar*df3.dyG
dnm
aa=dnm.cumsum('YC')

