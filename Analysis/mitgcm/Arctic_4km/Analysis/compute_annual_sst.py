import numpy as np
import os

root_folder = '/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/'
                                                                       
# expid = 'Exp02_0';                                                     
# expid = 'Exp02_1';                                                     
expid = 'Exp02_2';                                                     
                                                                       
fyear = 1992                                                         
lyear = 2018                                                           
datadir = root_folder+expid                                            
# os.chdir(datadir)                                                    
prename = '3DArcticOcean_monthly_THETA_' 

fname = datadir + '/' + prename +'*.nc'
list=sorted(glob.glob(fname))

time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)

# df = xr.open_mfdataset(list)
df = xr.open_mfdataset(list,combine='by_coords')     
# df = xr.open_mfdataset(list,combine='by_coords',chunks={'j':384,'time':12,'i_g':420})     

sst = df.THETA[:,0,:,:]
sst = sst.groupby('time.year').mean('time')
ds = sst.to_dataset(name='sst')

fname = root_folder + expid + '_SST.nc'
ds.to_netcdf(fname)