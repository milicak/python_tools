import gsw

root_folder = '/okyanus/users/milicak/dataset/SOSE/'

project_name = '1_over_3_degree'

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
fnames = root_folder+project_name+'/'+'*Vvel*'
list = sorted(glob.glob(fnames))
dfv = xr.open_mfdataset(list)

fnames = root_folder+project_name+'/'+'*Theta*'
list = sorted(glob.glob(fnames))
dft = xr.open_mfdataset(list)

fnames = root_folder+project_name+'/'+'*Salt*'
list = sorted(glob.glob(fnames))
dfs = xr.open_mfdataset(list)

voltrV = df.VVEL*df.hFacS*df.dxG*df.drF

Trx = voltrV.sum(('XC'))
Trxsum = Trx.cumsum('Z')
Trxsummean = Trxsum.mean('time')

plt.pcolormesh(df2.YG,df2.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()
