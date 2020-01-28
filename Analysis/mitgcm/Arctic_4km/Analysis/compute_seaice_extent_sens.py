import numpy as np        
import numpy.ma as ma     
import glob               
import xarray as xr       
import sys                
import pandas as pd       
import os                 
from datetime import date 

root_folder = '/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/'
                                                                       
expid = 'Exp02_0';                                                     
# expid = 'Exp02_1';                                                     
# expid = 'Exp02_2';                                                     
                                                                       
fyear = 1992                                                           
# fyear = 1992                                                         
lyear = 2008                                                           
# lyear = fyear+1                                                      
datadir = root_folder+expid                                            
# os.chdir(datadir)                                                    
prename = '2DArcticOcean_monthly_SIarea_'                                             

# SI_extent = np.zeros((lyear-fyear)*12)                                      
SI_extent = []
for year in range(fyear,lyear):                                             
    fname = datadir+'/'+prename+np.str(year)+'_1-12.nc'        
    print(fname)                                                        
    ds1 = xr.open_dataset(fname, chunks={'i':500, 'j':500})['SIarea']   
    ds1.data[ds1.data>0.15] = 1                                         
    ds1.data[ds1.data<=0.15] = 0                                        
    si = ds1*ds1.rA                                                     
    SI_extent=np.append(SI_extent,si.sum(dim=['i','j']).compute()) 


SI_extent = np.reshape(SI_extent,[lyear-fyear,12])  
df = pd.DataFrame(SI_extent)                        
if expid=='Exp02_1':
    df.to_csv("SI_extent_Atlantic_warm_1992_2007")               
elif expid=='Exp02_2':
    df.to_csv("SI_extent_Pacific_warm_1992_2007")               
elif expid=='Exp02_0':
    df.to_csv("SI_extent_ctrl_1992_2007")               

