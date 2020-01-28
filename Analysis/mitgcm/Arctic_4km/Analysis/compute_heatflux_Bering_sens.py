import numpy as np        
import numpy.ma as ma     
import glob               
import xarray as xr       
import sys                
import pandas as pd       
import os                 
from datetime import date 

root_folder = '/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/'
                                                                       
# expid = 'Exp02_0';                                                     
# expid = 'Exp02_1';                                                     
expid = 'Exp02_2';                                                     
                                                                       
fyear = 1992                                                           
# fyear = 1992                                                         
lyear = 2008                                                           
# lyear = fyear+1                                                      
datadir = root_folder+expid                                            
# os.chdir(datadir)                                                    
prename = '3DArcticOcean_monthly_UTHMASS_' 

# heat transport at bering strait Cm^3/s
HT_bering = []
for year in range(fyear,lyear):                                             
    fname = datadir+'/'+prename+np.str(year)+'_1-12.nc'        
    print(fname)                                                        
    df = xr.open_dataset(fname, chunks={'i':500, 'j':500})   
    uht = np.copy(df.UTHMASS[:,:,659:688,1449]*df.hFacW[:,659:688,1449]
                  *df.drF*np.copy(df.dyC[659:688,1449]))
    HT_bering = np.append(HT_bering,uht.sum(axis=(1,2)))


HT_bering= np.reshape(HT_bering,[lyear-fyear,12])  
df = pd.DataFrame(HT_bering)                        
if expid=='Exp02_1':
    df.to_csv("HT_bering_Atlantic_warm_1992_2007")               
elif expid=='Exp02_2':
    df.to_csv("HT_bering_Pacific_warm_1992_2007")               
elif expid=='Exp02_0':
    df.to_csv("HT_bering_ctrl_1992_2007")               

