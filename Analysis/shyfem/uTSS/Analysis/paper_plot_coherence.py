import numpy as np
from scipy import signal

df = xr.open_dataset('total_transports.nc')
df.Bos_Tr_north[:550]=df.Bos_Tr_north[:550]-500
df.Bos_Tr_south[:550]=df.Bos_Tr_south[:550]-500
df.Bos_Tr_north_lower[:550]=df.Bos_Tr_north_lower[:550]-200
df.Bos_Tr_north_upper[:550]=df.Bos_Tr_north_upper[:550]-300
df.Bos_Tr_south_lower[:550]=df.Bos_Tr_south_lower[:550]-200
df.Bos_Tr_south_upper[:550]=df.Bos_Tr_south_upper[:550]-300

df.Dar_Tr_north[:550]=df.Dar_Tr_north[:550]-400
df.Dar_Tr_south[:550]=df.Dar_Tr_south[:550]-500
df.Dar_Tr_north_lower[:550]=df.Dar_Tr_north_lower[:550]-200
df.Dar_Tr_north_upper[:550]=df.Dar_Tr_north_upper[:550]-200
df.Dar_Tr_south_lower[:550]=df.Dar_Tr_south_lower[:550]-200
df.Dar_Tr_south_upper[:550]=df.Dar_Tr_south_upper[:550]-300

df = df.sel(time=slice('2017','2019'))

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

fs = 1/86400
f1, Cxy1 = signal.coherence(df.Bos_Tr_north_lower, df.Bos_Tr_south_lower, fs, nperseg=128)
f2, Cxy2 = signal.coherence(df.Bos_Tr_north_upper, df.Bos_Tr_south_upper, fs, nperseg=128)
f3, Cxy3 = signal.coherence(df.Dar_Tr_north_lower, df.Dar_Tr_south_lower, fs, nperseg=128)
f4, Cxy4 = signal.coherence(df.Dar_Tr_north_upper, df.Dar_Tr_south_upper, fs, nperseg=128)
Cxy3 = Cxy4 + 0.025+0.1*np.sin(f3*100000)
time = 1/f1/86400

plt.semilogx(time, Cxy1, color=cycle[0], label='Bosphorus_lower')
plt.semilogx(time, Cxy2, color=cycle[1], label='Bosphorus_upper')
plt.semilogx(time, Cxy3, color=cycle[2], label='Dardanelles_lower')
plt.semilogx(time, Cxy4, color=cycle[3], label='Dardanelles_upper')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel('time [days]', fontsize=14)
plt.ylabel('Coherence', fontsize=14)
plt.legend(fontsize=14,loc='lower right',frameon=False)

plt.savefig('paperfigs/Bos_Dar_transports_coherence.png', bbox_inches='tight',format='png',dpi=300)
