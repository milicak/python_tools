
root_folder = '/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/'

project_name = 'Exp01_0'

# fnames = project_name+root_folder+'UVELMASS_2007_1-12_01cyc.nc'
fnames = root_folder+project_name+'/'+'*VVELMASS*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)

voltrV = df.VVELMASS*df.dxG*df.drF

Trx = voltrV.sum(('XC'))
Trxsum = Trx.cumsum('Z')
Trxsummean = Trxsum.mean('time')

plt.pcolormesh(df.YG,df.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()

