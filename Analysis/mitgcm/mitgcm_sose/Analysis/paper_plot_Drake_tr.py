
Drake_Donohue = 173.3
Drake_Cunningham = 134
Drake_Donohue = np.tile(Drake_Donohue,(84,1));
Drake_Cunningham = np.tile(Drake_Cunningham,(84,1));

fname = 'Exp01_0s_Drake_Tr.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp02_0s_Drake_Tr.nc'
df2 = xr.open_dataset(fname)
fname = 'Exp03_0s_Drake_Tr.nc'
df3 = xr.open_dataset(fname)
fname = 'Exp04_0s_Drake_Tr.nc'
df4 = xr.open_dataset(fname)
fname = 'Exp05_0s_Drake_Tr.nc'
df5 = xr.open_dataset(fname)

fig, ax = plt.subplots(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
# ax1 = plt.subplot(2, 1, 1)
ax1.plot(df.time,df1.DrakeTr*1e-6, 'b', label='Ctrl')
ax1.plot(df.time,df2.DrakeTr*1e-6, 'r', label='Wind')
ax1.plot(df.time,df3.DrakeTr*1e-6, color='orange', label='WindThermal')
ax1.plot(df.time,df4.DrakeTr*1e-6, 'g', label='BCCtrl')
ax1.plot(df.time,df5.DrakeTr*1e-6, 'm', label='BCWindThermal')
ax1.plot(df.time,Drake_Donohue, 'k--', label='Donohue et al. (2016)')
ax1.plot(df.time,Drake_Cunningham, 'k-.', label='Cunningham et al. (2003)')
# ax1.legend(fontsize=16)
ax1.legend()
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.set_ylabel('Transport [Sv]', fontsize=16)
ax1.set_xlabel('years', fontsize=16)

plt.savefig('paperfigs/Drake_Transport.png', bbox_inches='tight',format='png',dpi=300)

# ax2 = plt.subplot(2, 1, 2)
# ax2.plot(df.time,(df2.DrakeTr-df1.DrakeTr)*1e-6, 'r', label='Wind')
# ax2.plot(df.time,(df3.DrakeTr-df1.DrakeTr)*1e-6, color='orange', label='WindThermal')
# ax2.plot(df.time,(df4.DrakeTr-df1.DrakeTr)*1e-6, 'g', label='BCCtrl')
# ax2.plot(df.time,(df5.DrakeTr-df1.DrakeTr)*1e-6, 'k', label='BCWindThermal')
# ax2.legend(fontsize=16)
# ax2.tick_params(axis='x', labelsize=16)
# ax2.tick_params(axis='y', labelsize=16)
# ax2.set_ylabel('Transport [Sv]', fontsize=16)
# ax2.set_xlabel('years', fontsize=16)
