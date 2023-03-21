import numpy as np

df1 = xr.open_dataset('Bering_heat_transport.nc')
df2 = xr.open_dataset('Fram_heat_transport.nc')
df3 = xr.open_dataset('Davis_heat_transport.nc')
df4 = xr.open_dataset('Barents_heat_transport.nc')
df5 = xr.open_dataset('Fram_AW_heat_transport.nc')


total_heat_tr = df1.Bering_heat_tr-df2.Fram_heat_tr-df3.Davis_heat_tr+df4.Barents_heat_tr


# Positive values are INTO the Arctic Ocean
fig1, f1_axes = plt.subplots(figsize=(8,6),ncols=1, nrows=1,
                             constrained_layout=True, squeeze=False)
f1_axes[0,0].plot(df1['Time'],df1.Bering_heat_tr,label='Bering')
f1_axes[0,0].plot(df1['Time'],-df2.Fram_heat_tr,label='Fram')
f1_axes[0,0].plot(df1['Time'],-df3.Davis_heat_tr,label='Davis')
f1_axes[0,0].plot(df1['Time'],df4.Barents_heat_tr,label='Barents')
f1_axes[0,0].plot(df1['Time'],df5.Fram_AW_heat_tr,'k',label='Fram_AW')
f1_axes[0,0].set_ylim(-10,140);
f1_axes[0,0].set_xticklabels(['1995', "1996", "2000", "2004",'2008','2012','2016'], fontsize=14)
f1_axes[0,0].set_yticklabels(["-1","0", "20", "40", "60",'80','100','120','140','160'], fontsize=14)
f1_axes[0,0].set_ylabel(r'Heat transport [TW',fontsize = 14.0)
f1_axes[0,0].set_xlabel(r'Year',fontsize = 14.0)
f1_axes[0,0].grid()
f1_axes[0,0].legend()

plt.savefig('paperfigs/Heat_transports.png', bbox_inches='tight',format='png',dpi=300)
