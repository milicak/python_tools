import numpy as np


df = xr.open_dataset('Florida_Strait_transport.nc')
df_annual = df.groupby("time.year").mean("time")

# Positive values are INTO the Arctic Ocean
fig1, f1_axes = plt.subplots(figsize=(8,6),ncols=1, nrows=1,
                             constrained_layout=True, squeeze=False)
f1_axes[0,0].plot(df['time'],df.Florida_tr,label='monthly')
f1_axes[0,0].plot(df['time'][5::12],df_annual.Florida_tr,label='annual')
f1_axes[0,0].set_ylim(10,40);
f1_axes[0,0].grid()
f1_axes[0,0].set_xticklabels(['1995', "1996", "2000", "2004",'2008','2012','2016'], fontsize=14)
f1_axes[0,0].set_yticklabels(["10","15", "20", "25", "30",'35','40'], fontsize=14)
f1_axes[0,0].set_ylabel(r'Florida Strait volume transport [Sv]',fontsize = 14.0)
f1_axes[0,0].set_xlabel(r'Year',fontsize = 14.0)
f1_axes[0,0].legend()
f1_axes[0,0].text(df.time[0],38,'mean = 24.13 Sv',fontsize = 14.0)

plt.savefig('paperfigs/Florida_Volume_transports.png', bbox_inches='tight',format='png',dpi=300)

