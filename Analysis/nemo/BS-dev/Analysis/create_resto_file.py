import numpy as np
import scipy.io                                            

df = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')

nk = df.z.shape[0]
ni = df.x.shape[0]
nj = df.y.shape[0]

# mask marmara sea
tau1 = 1./(1*86400) 
tau2 = 1./(5*86400) 

xind1 = 40
xind2 = 85
yind1 = 5
yind2 = 22

tau = np.zeros((nk, nj, ni))
Lend = yind2-yind1
taunew = np.zeros((Lend))
for ind in range(0,yind2-yind1):
    taunew[ind] = tau1-(tau1-tau2)*(ind+1)/Lend

tnew = np.tile(taunew[None, :, None], (121, 1, xind2-xind1))
tau[:,yind1:yind2,xind1:xind2] = tnew

fout  = 'resto.nc'
rg = scipy.io.netcdf_file(fout,'w')                            
# Dimensions                                                               
rg.createDimension('z',nk)                                               
rg.createDimension('x',ni)                                               
rg.createDimension('y',nj)                                               
# Variables                                                                
tempvar  = rg.createVariable('resto','float32',('z','y','x',)) 
# Values                                                                   
tempvar[:] = tau                                         
rg.close()                                                                 
