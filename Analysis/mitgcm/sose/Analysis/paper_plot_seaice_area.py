import numpy as np


# load the data
df1=xr.open_dataset('Exp01_0_SIarea.nc')
df2=xr.open_dataset('Exp02_0_SIarea.nc')
df3=xr.open_dataset('Exp03_0_SIarea.nc')

d1 = df1.SIarea
d2 = df2.SIarea
d3 = df3.SIarea

# time = pd.date_range("2005-01-05", freq="5D", periods=511)
# d1['time'] = time
# d2['time'] = time
# d3['time'] = time
# d1 = d1.resample(time="1MS").mean('time')
# d2 = d2.resample(time="1MS").mean('time')
# d3 = d3.resample(time="1MS").mean('time')
#
fig, ax = plt.subplots(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
# ax1 = plt.subplot(2, 1, 1)
ax1.plot(d1.time,d1*1e-12, 'b', label='Ctrl')
ax1.plot(d2.time,d2*1e-12, 'r', label='Wind')
ax1.plot(d3.time,d3*1e-12, color='orange', label='WindThermal')
ax1.tick_params(labelsize=16)
ax1.set_ylabel('Sea Ice Area [10$^6$ km$^2$]', fontsize=16)
ax1.set_xlabel('years', fontsize=16)
plt.grid()
ax1.legend(loc='lower right',frameon=False,bbox_to_anchor=(1.23, 0.85));

plt.savefig('paperfigs/seaice_area.png', bbox_inches='tight',format='png',dpi=300)
