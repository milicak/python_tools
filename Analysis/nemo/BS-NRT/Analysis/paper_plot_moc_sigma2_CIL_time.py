import numpy as np
import os
# read temperature change
# ls1=sorted(glob.glob('/work/opa/da01720/Tools/hov/rean16/*domain.nc'))
ls1=sorted(glob.glob('/work/opa/da01720/Tools/hov/rean16_2020/*domain.nc'))
df=xr.open_mfdataset(ls1)   
# time = pd.date_range("1993-01-01", freq="D", periods=9861)  
time = pd.date_range("1993-01-01", freq="D", periods=10227)  
df['time']=time 
dfm=df.resample(time="1MS").mean(dim="time")
# 70 meter depth which k = 12
dfm70=dfm.votemper[:,12]                                                                                                         
dfm70y=dfm70.groupby('time.year').mean('time')    

# total volume
base='/work/opa/sc33616/nemo/tools/OC/daily_rean16_vol/'
years=np.arange(1993,2021,1)

buffer=[]
for y in years:
    buffer.append(np.nanmean(np.load(os.path.join(base,f"vol_{y}.npy"))))

buffer=np.array(buffer) 

# plotyy
# df2=xr.open_dataset('maxmocy.nc') 
df2=xr.open_dataset('maxmocy_1993_2020.nc') 

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, ax1 = plt.subplots()
fig.subplots_adjust(right=0.75)
# color = 'tab:red'
color = cycle[3]
ax1.set_ylabel(r'MOC in $\sigma_2$ space [Sv]', color=color, fontsize=14) 
ax1.plot(df2.year, df2.vol_sigma_tr, '-*', color=color)
ax1.set_xlabel('year', fontsize=14)
ax1.tick_params(axis='y', labelcolor=color)
# plt.gca().invert_yaxis()  
ax1.set_yticklabels(ax1.get_yticks(), rotation=0, fontsize=14);
ax1.set_yticklabels(['0.07', '0.08', '0.09', '0.1', '0.11', '0.12', '0.13'], fontsize=14); 
yticks = np.arange(0.07, 0.132, 0.01)  
ax1.set_yticks(yticks);
ax1.set_xticklabels(['1995', '2000', '2005', '2010', '2015', '2020'], fontsize=14); 
xticks = np.arange(1995, 2025, 5)  
ax1.set_xticks(xticks);


ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

df3 = xr.open_dataset('/data/opa/bs-mod/osr6/Black_Sea_CIL_Cold_Content_Annual.nc')   
df3 = df3.sel(time=slice('1993','2019'))  
ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax3.spines["right"].set_position(("axes", 1.2))

color = cycle[0]
# color = 'tab:blue'
# ax2.set_ylabel(r'Temperature [$^\circ$C] at 70m', color=color, fontsize=14) 
ax2.set_ylabel(r'Volume of Temperature [$<8.35^\circ$C] [1000xkm^3]', color=color, fontsize=14) 
# ax2.plot(np.int16(df2.year), dfm70y,'-*' , color=color)
ax2.plot(np.int16(df2.year), buffer*1e-12,'-*' , color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_yticklabels(ax2.get_yticks(), rotation=0, fontsize=14);

color = cycle[1]
ax3.set_ylabel(r'CIL cold content [M J/m$^2$]', color=color, fontsize=14) 
ax3.plot(df2.year[:-1], np.copy(df3.C_Argo), '-', color=color) 
ax3.plot(df2.year[:-1], np.copy(df3.C_ShipCast), '-', color=color) 
ax3.tick_params(axis='y', labelcolor=color)
ax3.set_yticklabels(ax3.get_yticks(), rotation=0, fontsize=14);


# plt.savefig('paperfigs/MOC_sigma2_vs_CIL_vs_70mtemp.png', bbox_inches='tight',format='png',dpi=300)
plt.savefig('paperfigs/MOC_sigma2_vs_CIL_vs_voltemp.jpeg',
            bbox_inches='tight',format='jpeg',dpi=300)
# plt.savefig('paperfigs/MOC_sigma2_vs_CIL_vs_voltemp.png', bbox_inches='tight',format='png',dpi=300)
