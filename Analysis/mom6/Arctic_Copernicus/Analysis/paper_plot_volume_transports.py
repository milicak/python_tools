import numpy as np

df1 = xr.open_dataset('Bering_volume_transport.nc')
df2 = xr.open_dataset('Fram_volume_transport.nc')
df3 = xr.open_dataset('Davis_volume_transport.nc')
df4 = xr.open_dataset('Barents_volume_transport.nc')


total_volume_tr = df1.Bering_volume_tr-df2.Fram_volume_tr-df3.Davis_volume_tr+df4.Barents_V_volume_tr


# Positive values are INTO the Arctic Ocean
fig1, f1_axes = plt.subplots(figsize=(8,6),ncols=1, nrows=1,
                             constrained_layout=True, squeeze=False)
f1_axes[0,0].plot(df1['Time'],df1.Bering_volume_tr,label='Bering')
f1_axes[0,0].plot(df1['Time'],-df2.Fram_volume_tr,label='Fram')
f1_axes[0,0].plot(df1['Time'],-df3.Davis_volume_tr,label='Davis')
f1_axes[0,0].plot(df1['Time'],df4.Barents_V_volume_tr,label='Barents')
# f1_axes[0,0].plot(df1['Time'],total_volume_tr,'k',label='Barents')
f1_axes[0,0].set_ylim(-6,6);
f1_axes[0,0].set_xticklabels(['1995', "1996", "2000", "2004",'2008','2012','2016'], fontsize=14)
f1_axes[0,0].set_yticklabels(["-6","-4", "-2", "0", "2",'4','6'], fontsize=14)
f1_axes[0,0].set_ylabel(r'Volume transport [Sv]',fontsize = 14.0)
f1_axes[0,0].set_xlabel(r'Year',fontsize = 14.0)
f1_axes[0,0].grid()
f1_axes[0,0].legend()

plt.savefig('paperfigs/Volume_transports.png', bbox_inches='tight',format='png',dpi=300)
