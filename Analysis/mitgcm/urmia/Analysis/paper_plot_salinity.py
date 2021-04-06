import xmitgcm

time = 100

df = xmitgcm.open_mdsdataset('/okyanus/users/milicak/models/MITgcm/Projects/urmia/work/',
                             prefix='SALT',geometry='cartesian',read_grid='True',
                             ref_date='2018-12-31 12:0:0',delta_t=2)


df2 = xmitgcm.open_mdsdataset('/okyanus/users/milicak/models/MITgcm/Projects/urmia/work/',
                        prefix='S',geometry='cartesian',read_grid='True')

mask = np.double(S1.maskC)
kmax = np.int16(np.sum(mask,0))
kmax = xr.DataArray(kmax, dims={'points1','points2'})

Sbottom = np.zeros((df2.S.shape[2],df2.S.shape[3]))
for ii in range(0,df2.S.shape[3]):
    for jj in range(0,df2.S.shape[2]):
        Sbottom[jj,ii] = df.SALT[time,kmax[jj,ii],jj,ii]
