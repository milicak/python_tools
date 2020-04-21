import numpy as np        
import numpy.ma as ma     
import glob               
import xarray as xr       
import sys                
import pandas as pd       
import os                 
from datetime import date 
import datetime

root_folder = '/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/'
                                                                       
# expid = 'Exp02_0';                                                     
# expid = 'Exp02_1';                                                     
expid = 'Exp02_2';                                                     
                                                                       
fyear = 1992                                                           
# fyear = 1992                                                         
lyear = 2018                                                           
# lyear = fyear+1                                                      
datadir = root_folder+expid                                            
# os.chdir(datadir)                                                    
prename = '3DArcticOcean_monthly_UTHMASS_' 

fname = datadir + '/' + prename +'*.nc'
list=sorted(glob.glob(fname))

time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)

# df = xr.open_mfdataset(list)
df = xr.open_mfdataset(list,combine='by_coords')     
# old way 
time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)
# new way
# arr = np.array([datetime.datetime(1991, 12, 16) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = np.array([datetime.datetime(1992, 1, 1) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = arr[:300]
newtime = xr.Dataset({'time': arr}) 
# df['time'] = newtime['time']

# heat transport at fram strait Cm^3/s
ds = df.UTHMASS*df.dyG*df.drF

ht = ds.data[:,:,496:664,580]
HT_fram = ht.sum(axis=(1,2))
dnm = HT_fram.compute()

# data = xr.DataArray(dnm, dims=('time'), coords={'time': time})
data = xr.DataArray(dnm, dims=('time'), coords={'time': df['time']})
dsnew = data.to_dataset(name='heat_transport')
fname = expid + '_Fram_heat_transport.nc'
dsnew.to_netcdf(fname)

# HT_fram = []
# for year in range(fyear,lyear):                                             
#     fname = datadir+'/'+prename+np.str(year)+'_1-12.nc'        
#     print(fname)                                                        
#     df = xr.open_dataset(fname, chunks={'i':500, 'j':500})   
#     uht = np.copy(df.UTHMASS[:,:,496:664,580]*df.hFacW[:,496:664,580]
#                   *df.drF*np.copy(df.dyC[496:664,580]))
#     HT_fram = np.append(HT_fram,uht.sum(axis=(1,2)))
#
#
# HT_fram= np.reshape(HT_fram,[lyear-fyear,12])  
# df = pd.DataFrame(HT_fram)                        
# if expid=='Exp02_1':
#     df.to_csv("HT_fram_Atlantic_warm_1992_2007")               
# elif expid=='Exp02_2':
#     df.to_csv("HT_fram_Pacific_warm_1992_2007")               
# elif expid=='Exp02_0':
#     df.to_csv("HT_fram_ctrl_1992_2007")               
#
