import numpy as np

Drake_Donohue = 173.3
Drake_Cunningham = 134
Drake_Donohue = np.tile(Drake_Donohue,(85,1));
Drake_Cunningham = np.tile(Drake_Cunningham,(85,1));

# load the data
df1=xr.open_dataset('Exp01_0_Drake_Tr.nc')
df2=xr.open_dataset('Exp02_0_Drake_Tr.nc')
df3=xr.open_dataset('Exp03_0_Drake_Tr.nc')

d1 = df1.DrakeTr[1:]
d2 = df2.DrakeTr[1:]
d3 = df3.DrakeTr[1:]
time = pd.date_range("2005-01-05", freq="5D", periods=511)
d1['time'] = time
d2['time'] = time
d3['time'] = time
d1 = d1.resample(time="1MS").mean('time')
d2 = d2.resample(time="1MS").mean('time')
d3 = d3.resample(time="1MS").mean('time')

fig, ax = plt.subplots(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
# ax1 = plt.subplot(2, 1, 1)
ax1.plot(d1.time,d1*1e-6, 'b', label='Ctrl')
ax1.plot(d1.time,d2*1e-6, 'r', label='Wind')
ax1.plot(d1.time,d3*1e-6, color='orange', label='WindThermal')
ax1.plot(d1.time,Drake_Donohue[1:], 'k--', label='Donohue et al. (2016)')
ax1.plot(d1.time,Drake_Cunningham[1:], 'k-.', label='Cunningham et al. (2003)')
# ax1.legend(fontsize=16)
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.set_ylabel('Transport [Sv]', fontsize=16)
ax1.set_xlabel('years', fontsize=16)
plt.grid()
ax1.set_ylim(80,200);
ax1.legend(loc='lower right',frameon=False);

plt.savefig('paperfigs/Drake_Transport_monthly.png', bbox_inches='tight',format='png',dpi=300)
