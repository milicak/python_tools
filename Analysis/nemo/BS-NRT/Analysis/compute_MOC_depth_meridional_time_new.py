import numpy as np


ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/nemo/BS/MOC_data/*meridional*sigma2*BSe3r1*')) 
ls1 = ls1[:-2]   

# from 1993 - 2020
time = pd.date_range("1993-01-01", freq="D", periods=10227)  
df = xr.open_dataset(ls1[0])

vol_sigma_tr = np.zeros((time.shape[0],df.lat.shape[0],df.sigma2_bin.shape[0]))

for ii in enumerate(ls1): 
    print(ii[0])
    df = xr.open_dataset(ii[1])
    vol_sigma_tr[ii[0],:,:] = np.copy(df.vol_sigma_tr)  


# create dataset
ds = xr.Dataset({'vol_sigma_tr': (('time','lat','sigma2_bin'),vol_sigma_tr)},
                {'time': time, 'lat': df.lat, 'sigma2_bin': df.sigma2_bin})

dsm = ds.resample(time="1MS").mean(dim="time")
