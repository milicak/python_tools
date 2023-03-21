import numpy as np

# ds1 = xr.open_dataset('Dardanelles_transports.nc')
# ds2 = xr.open_dataset('Bosphorus_transports.nc')
# time = pd.date_range('2016-01-01', freq='D', periods=1461)
#
# df = xr.Dataset(
#     {
#         "Dar_Tr_south": (("time"), ds1.Dar_Tr_south),
#         "Dar_Tr_south_lower": (("time"), ds1.Dar_Tr_south_lower),
#         "Dar_Tr_south_upper": (("time"), ds1.Dar_Tr_south_upper),
#         "Dar_Tr_north": (("time"), ds1.Dar_Tr_north),
#         "Dar_Tr_north_lower": (("time"), ds1.Dar_Tr_north_lower),
#         "Dar_Tr_north_upper": (("time"), ds1.Dar_Tr_north_upper),
#         "Bos_Tr_south": (("time"), ds2.Bos_Tr_south),
#         "Bos_Tr_south_lower": (("time"), ds2.Bos_Tr_south_lower),
#         "Bos_Tr_south_upper": (("time"), ds2.Bos_Tr_south_upper),
#         "Bos_Tr_north": (("time"), ds2.Bos_Tr_north),
#         "Bos_Tr_north_lower": (("time"), ds2.Bos_Tr_north_lower),
#         "Bos_Tr_north_upper": (("time"), ds2.Bos_Tr_north_upper),
#     },
#     {"time": time},
# )
#

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
ds = df.resample(time='1MS').mean('time')
ds['time'] = ds.time + pd.Timedelta("15 day")

fig1, f1_axes = plt.subplots(figsize=(8,12),ncols=1, nrows=4, constrained_layout=True)
f1_axes[0].plot(df.time,0.8*df.Bos_Tr_north,color='lightgrey')
f1_axes[0].plot(df.time,df.Bos_Tr_north_upper,color='lightskyblue')
f1_axes[0].plot(df.time,0.8*df.Bos_Tr_north_lower,color='lightcoral')
f1_axes[0].plot(ds.time,ds.Bos_Tr_north_upper,color='b',label='UL')
f1_axes[0].plot(ds.time,0.8*ds.Bos_Tr_north_lower,color='r',label='LL')
f1_axes[0].plot(ds.time,0.8*ds.Bos_Tr_north,color='k',label='NET')
f1_axes[0].set_ylim(-2000,2000);
f1_axes[0].grid()
f1_axes[0].set_ylabel(r'Transport [km$^3$ yr $^{-1}$]')
f1_axes[0].legend(fontsize=14,loc='lower left',frameon=False,ncol=3)
# f1_axes[0].legend(fontsize=14,loc='upper right',frameon=False)
# f1_axes[0].tick_params(labelsize=14)
f1_axes[0].text(ds.time[0]-pd.Timedelta("60 day"),-1900,'a)',fontsize=14);
# fig2
f1_axes[1].plot(df.time,0.8*df.Bos_Tr_south,color='lightgrey')
f1_axes[1].plot(df.time,df.Bos_Tr_south_upper,color='lightskyblue')
f1_axes[1].plot(df.time,0.8*df.Bos_Tr_south_lower,color='lightcoral')
f1_axes[1].plot(ds.time,0.8*ds.Bos_Tr_south,color='k')
f1_axes[1].plot(ds.time,ds.Bos_Tr_south_upper,color='b')
f1_axes[1].plot(ds.time,0.8*ds.Bos_Tr_south_lower,color='r')
f1_axes[1].set_ylim(-2000,2000);
f1_axes[1].grid()
f1_axes[1].set_ylabel(r'Transport [km$^3$ yr $^{-1}$]')
f1_axes[1].text(ds.time[0]-pd.Timedelta("60 day"),-1900,'b)',fontsize=14);

f1_axes[2].plot(df.time,0.8*df.Dar_Tr_north,color='lightgrey')
f1_axes[2].plot(df.time,0.8*df.Dar_Tr_north_upper,color='lightskyblue')
f1_axes[2].plot(df.time,0.8*df.Dar_Tr_north_lower,color='lightcoral')
f1_axes[2].plot(ds.time,0.8*ds.Dar_Tr_north_upper,color='b',label='UL')
f1_axes[2].plot(ds.time,0.8*ds.Dar_Tr_north_lower,color='r',label='LL')
f1_axes[2].plot(ds.time,0.8*ds.Dar_Tr_north,color='k',label='NET')
f1_axes[2].set_ylim(-2000,2000);
f1_axes[2].grid()
f1_axes[2].set_ylabel(r'Transport [km$^3$ yr $^{-1}$]')
f1_axes[2].legend(fontsize=14,loc='lower left',frameon=False,ncol=3)
f1_axes[2].text(ds.time[0]-pd.Timedelta("60 day"),-1900,'c)',fontsize=14);

f1_axes[3].plot(df.time,0.8*df.Dar_Tr_south,color='lightgrey')
f1_axes[3].plot(df.time,0.8*df.Dar_Tr_south_upper,color='lightskyblue')
f1_axes[3].plot(df.time,0.8*df.Dar_Tr_south_lower,color='lightcoral')
f1_axes[3].plot(ds.time,0.8*ds.Dar_Tr_south,color='k')
f1_axes[3].plot(ds.time,0.8*ds.Dar_Tr_south_upper,color='b')
f1_axes[3].plot(ds.time,0.8*ds.Dar_Tr_south_lower,color='r')
f1_axes[3].set_ylim(-2000,2000);
f1_axes[3].grid()
f1_axes[3].set_ylabel(r'Transport [km$^3$ yr $^{-1}$]')
f1_axes[3].text(ds.time[0]-pd.Timedelta("60 day"),-1900,'d)',fontsize=14);

# plt.savefig('paperfigs/Bos_Dar_transports.png', bbox_inches='tight',format='png',dpi=300)





