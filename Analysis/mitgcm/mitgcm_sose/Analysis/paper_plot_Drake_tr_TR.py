
Drake_Donohue = 173.3
Drake_Cunningham = 134
Drake_Donohue = np.tile(Drake_Donohue,(84,1));
Drake_Cunningham = np.tile(Drake_Cunningham,(84,1));

fname = 'Exp01_0s_Drake_Tr.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp03_0s_Drake_Tr.nc'
df3 = xr.open_dataset(fname)

fig, ax = plt.subplots(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
# ax1 = plt.subplot(2, 1, 1)
ax1.plot(df1.time,df1.DrakeTr*1e-6, 'b', label='Kontrol')
ax1.plot(df1.time,df3.DrakeTr*1e-6, 'r', label=r'1.5$\times$rüzgar')
ax1.plot(df1.time,Drake_Donohue, 'k--', label='Donohue vd. (2016)')
ax1.plot(df1.time,Drake_Cunningham, 'k-.', label='Cunningham vd. (2003)')
ax1.legend(fontsize=16)
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.set_ylabel('Taşınım [Sv]', fontsize=16)
ax1.set_xlabel('Yıllar', fontsize=16)

plt.savefig('paperfigs_tr/Drake_Transport.png', bbox_inches='tight',format='png',dpi=300)

