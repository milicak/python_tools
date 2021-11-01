import numpy as np



ls1 = sorted(glob.glob('/work/opa/sc02915/exp/bs-nrt_8.9/rebuilt/*2014*grid_T*'))
df1 = xr.open_mfdataset(ls1)   
d1 = df1.so.mean('time_counter') 

ls1=sorted(glob.glob('/data/opa/bs-mod/experiments/bs-nrt_8.14_nobbl/rebuilt/*grid_T*'))                                                      
df1 = xr.open_mfdataset(ls1)  
d1 = df1.thetao.mean('time_counter') 

ls1=sorted(glob.glob('/data/opa/bs-mod/experiments/bs-nrt_8.22/rebuilt/*grid_T*'))
df3 = xr.open_mfdataset(ls1)  
d3 = df3.thetao.mean('time_counter') 

y1 = np.array([20,21,21,22,22,23,24,25,26,27,28,28,29,30,31,32,33,34,35,36,37])                                                
x1 = np.array([70,71,71,71,72,72,72,72,72,72,73,74,74,74,74,74,74,74,74,74,74]) 
y1 = xr.DataArray(y1, dims='points')    
x1 = xr.DataArray(x1, dims='points')    



s1 = d1.isel(y=y1,x=x1)
s2 = d2.isel(y=y1,x=x1)
s3 = d3.isel(y=y1,x=x1)
