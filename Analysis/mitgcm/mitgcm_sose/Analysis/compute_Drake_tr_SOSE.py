


root_folder = '/okyanus/users/milicak/dataset/SOSE/'

project_name = '1_over_3_degree'

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
fnames = root_folder+project_name+'/'+'*Uvel*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)

# location of the Drake Passage for 1/3 degree SOSE simulation
iind = 883
jind = slice(119, 191, 1)

voltr = df.UVEL*df.hFacW*df.dyG*df.drF
Tr = voltr[:,:,jind,iind]
Trtime = Tr.sum(('Z','YC'))


