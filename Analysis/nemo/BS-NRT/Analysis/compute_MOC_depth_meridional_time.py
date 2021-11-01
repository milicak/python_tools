import numpy as np

ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/nemo/BS/MOC_data/*depth*meridional*BSe3r1*'))
df=xr.open_mfdataset(ls1)       
dfm = df.resample(time="1MS").mean(dim="time")  

dfy = df.groupby('time.year').mean(dim="time")
dfy_surf=dfy.MOC_meridional[:,:21,:]
MOC=dfy_surf.max(('lat','depth'))*1e-6   

# first 250 meter
dfm_surf = dfm.MOC_meridional[:,:21,:]  
MOC250m = dfm_surf.max(('lat','depth'))*1e-6   
MOC250m = MOC250m.groupby('time.year').mean('time') 

# first 70 meter
dfm_surf = dfm.MOC_meridional[:,:12,:]  
MOC70m = dfm_surf.max(('lat','depth'))*1e-6   
MOC70m = MOC70m.groupby('time.year').mean('time') 

# max MOC 
MOC = dfm.max(('lat','depth'))*1e-6   
MOC = MOC.groupby('time.year').mean('time') 

ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/nemo/BS/MOC_data/*meridional*sigma2*BSe3r1*')) 
df = xr.open_mfdataset(ls1,concat_dim='time') 
# from 1993 - 2019
# time = pd.date_range("1993-01-01", freq="D", periods=9861)  
# from 1993 - 2020
time = pd.date_range("1993-01-01", freq="D", periods=10227)  
df['time'] = time
# correction because of the Z error in the computation
df['vol_sigma_tr'] = -df.vol_sigma_tr
dfm = df.resample(time="1MS").mean(dim="time")
maxmoc = dfm.vol_sigma_tr[:,:,55:65].max(('lat','sigma2_bin'))*1e-6 
maxmocy = maxmoc.groupby('time.year').mean('time') 
maxmocy = maxmocy.to_dataset(name='vol_sigma_tr')
maxmocy.to_netcdf('maxmocy_1993_2020.nc')


# read temperature change
ls1=sorted(glob.glob('/work/opa/da01720/Tools/hov/rean16/*domain.nc'))
df=xr.open_mfdataset(ls1)   
time = pd.date_range("1993-01-01", freq="D", periods=9861)  
df['time']=time 
dfm=df.resample(time="1MS").mean(dim="time")
# 70 meter depth which k = 12
dfm70=dfm.votemper[:,12]                                                                                                         
dfm70y=dfm70.groupby('time.year').mean('time')    

# plotyy
df2=xr.open_dataset('maxmocy.nc') 

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('year')
ax1.set_ylabel('temperature at 70m', color=color)
ax1.plot(df2.year, dfm70y,'-*' , color=color)
ax1.tick_params(axis='y', labelcolor=color)
plt.gca().invert_yaxis()  

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_ylabel(r'MOC($\sigma_2$)', color=color) 
ax2.plot(df2.year, df2.vol_sigma_tr, '-*', color=color)
ax2.tick_params(axis='y', labelcolor=color)



# 10 level density is around 23.65-23.9

# aa=dfm.vol_sigma_tr[:,:,57:62].sum('sigma2_bin').max('lat')*1e-6  

dfy = df.groupby('time.year').mean(dim="time") 
dnm = np.copy(dfy.vol_sigma_tr) 
cc = dnm[:,:,50:65].sum(axis=2) 
plt.plot(np.arange(1993,2020),cc.min(axis=1)*1e-6,'-*')
