import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
# from mpl_toolkits.basemap import Basemap                                            
from scipy.interpolate import interp1d
import geopy.distance

plt.ion()

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

# expid = 'Exp01.2'
# expid = 'Exp_20160101'
# expid = 'Exp_2016_analysis'
# expid = 'Exp_2016_analysis_keps'
expid = 'Exp_2016_analysis_newTSIC'

df = pd.read_excel('20190220_1033_ODV_TMAMVerisi.xlsx',sheetname='2018data')

timecol = 'yyyy-mm-ddThh:mm:ss.sss'
depthcol = 'Depth [m]'
tempcol = 'Temperature [degC]'
saltcol = 'Salinity [psu]'

stations = df[timecol].unique()

# for each station
timeind = 0
# expname = '2018-01_BUTUNLESIK_MARMARA'
# expname = '2018-04_BUTUNLESIK_MARMARA'
expname = '2018-08_BUTUNLESIK_MARMARA'
tmp_pro = df.loc[df['Cruise']==expname]

# stationame = 'M23'  # 'M23' 'B2' 'B2B' 'B5' 'B7' 'K0A' 'K0' 'M14' 'M18' 'K0' 'K2' 'BL1' 'B14' 'B13' 'M1' 'M8' 'MBA' 'M20'
# stationames = ['M23', 'B2', 'B2B', 'B5', 'B7', 'K0A', 'K0', 'M14']
# stationames = ['M20', 'M23', 'M8', 'M14', 'MY2', 'MBA', 'M1', 'SB', 'B7',
               # 'B13', 'B14']

stationames = ['K0', '45C', 'MD102', 'M20', 'MD104', 'M14A', 'YSA', 'DIPTAR1Y',
               'DIPTAR2Y', 'DIPTAR3Y', 'DIPTAR4Y', 'MD75', 'AR1', 'MD18',
               'MD101',  'MD13A', 'KD1', 'MD67', 'MD10A', 'MD10B']
 
 


for stationame in stationames:
    print(stationame)
    tmp_prof = tmp_pro.loc[tmp_pro['Station']==stationame]
    # compute the day of the year
    year = tmp_prof[timecol].tolist()[timeind][0:4]
    mnth = tmp_prof[timecol].tolist()[timeind][5:7]
    day = tmp_prof[timecol].tolist()[timeind][8:10]
    date = pd.to_datetime(year+mnth+day, format='%Y%m%d')
    new_year_day = pd.Timestamp(year=2016, month=1, day=1)
    day_of_the_year = (date - new_year_day).days + 1
    ind_name = 'stations_indices/' + stationame + '_obs.csv'
    di = pd.read_csv(ind_name)
    d_utss = np.copy(di['d_utss'])
    inds_utss = np.copy(di['inds_utss'])
    
    
    sdate = "%4.4d" % (day_of_the_year)
    fname = root_folder+project_name+'/'+expid+'/OUT_2017_2018_2019/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
    ds = xr.open_dataset(fname)
    utss_temp = ds.temperature[0,:,:]
    utss_salt = ds.salinity[0,:,:]
    mask = utss_temp.where(utss_temp==0,1,0)

    L = 0.25 # correlation distance
    
    ### calculate weights from d_mfs, d_utss, d_bsfs
    w_utss = np.exp(-(d_utss/L)**2)
    
    ## I have to expand weights array to make
    ## multiplication. Basically I have to repeat the
    ## weights array along depth direction making this
    ## transformation of shape
    ## (nodes,neighbours) --> (nodes,neighbours,depth)
    
    w_utss = np.expand_dims(w_utss,axis=2)
    w_utss = np.repeat(w_utss,ds.level.shape[0],axis=1)
    
    # give a minimum value to weights
    w_utss = np.clip(w_utss,a_min=0.01,a_max=1)
    
    w_utss = np.ma.masked_where( utss_temp[inds_utss,:] == 0, w_utss)
    
    ######################################################
    ### CALCULATE WEIGHTED AVERAGE for T,Sprint 'Calculating weighted average..'
    utss_part_temp = np.sum( w_utss * utss_temp[inds_utss], axis = 0) / np.sum( w_utss, axis = 0)
    utss_part_salt = np.sum( w_utss * utss_salt[inds_utss], axis = 0) / np.sum( w_utss, axis = 0)

    
    df1 = pd.DataFrame({"Temp_obs": np.array(tmp_prof[tempcol]),"Salt_obs": np.array(tmp_prof[saltcol]),
                        "zr_obs": -np.array(tmp_prof[depthcol])})
    
    df2 = pd.DataFrame({"Temp_uTSS": np.copy(utss_part_temp),"Salt_uTSS":
                        np.copy(utss_part_salt),
                       "zr_uTSS": -ds.level.compute()})
    
    fname = 'obs_ctd_netcdf/2018_08_'+stationame+'_obs.csv'
    df1.to_csv(fname)
    fname = 'obs_ctd_netcdf/2018_08_'+stationame+'_uTSS.csv'
    df2.to_csv(fname)




