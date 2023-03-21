import xmitgcm

root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
outname = project_name + '_stream_function.nc'
fnames = root_folder+project_name+'/'

fnames = root_folder+project_name+'/*UVELMASS*'
grname = root_folder+project_name+'/grid.nc'
list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)
gr = xr.open_dataset(grname)

mask = gr.maskC[0,:,:]
mask = np.multiply(mask,1)
mask = mask.rename({"XC": "XG"})

df = df.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
Ubar = df.UVELMASS*gr.drF


Ubar = Ubar.sum('Z')
dnm = Ubar*gr.dyG
dnm
aa = dnm.cumsum('YC')
aa = aa.mean('time')
dnm = aa*mask.data
ds = dnm.to_dataset(name='stream_function')
ds.to_netcdf(outname)




