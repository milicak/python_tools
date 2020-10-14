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

df = pd.read_excel('20190220_1033_ODV_TMAMVerisi.xlsx',sheetname='2017data')

timecol = 'yyyy-mm-ddThh:mm:ss.sss'
depthcol = 'Depth [m]'
tempcol = 'Temperature [degC]'
saltcol = 'Salinity [psu]'

stations = df[timecol].unique()

# for each station
timeind = 0
expname = '2017-01_ISKI_MARMARA'
tmp_pro = df.loc[df['Cruise']==expname]
# stationame = 'M23'  # 'M23' 'B2' 'B2B' 'B5' 'B7' 'K0A' 'K0' 'M14' 'M18' 'K0' 'K2' 'BL1' 'B14' 'B13' 'M1' 'M8' 'MBA' 'M20'
# stationames = ['M23', 'B2', 'B2B', 'B5', 'B7', 'K0A', 'K0', 'M14']
stationames = ['M20', 'M23', 'M8', 'M14', 'MY2', 'MBA', 'M1', 'SB', 'B7',
               'B13', 'B14']
for stationame in stationames:
    tmp_prof = tmp_pro.loc[tmp_pro['Station']==stationame]
    # tmp_pro = df.loc[df['Cruise']=='2016-09_ISKI_MARMARA']
    # tmp_prof = tmp_pro.loc[tmp_pro['Station']=='BL1']
    # tmp_prof = tmp_pro.loc[tmp_pro['Station']=='K0']
    # tmp_prof = tmp_pro.loc[tmp_pro['Station']=='K2']
    # compute the day of the year
    year = tmp_prof[timecol].tolist()[timeind][0:4]
    mnth = tmp_prof[timecol].tolist()[timeind][5:7]
    day = tmp_prof[timecol].tolist()[timeind][8:10]
    date = pd.to_datetime(year+mnth+day, format='%Y%m%d')
    new_year_day = pd.Timestamp(year=2016, month=1, day=1)
    day_of_the_year = (date - new_year_day).days + 1
    day_of_the_year
    
    lon_thlweg = np.array([ tmp_prof['Longitude'].tolist()[0]])
    lat_thlweg = np.array([ tmp_prof['Latitude'].tolist()[0]])
    
    sdate = "%4.4d" % (day_of_the_year)
    fname = root_folder+project_name+'/'+expid+'/OUT_2017_2018_2019/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
    fname
    
    domask = False
    # create locstream 
    coord_sys = ESMF.CoordSys.SPH_DEG
    locstream = ESMF.LocStream(lon_thlweg.shape[0], name="uTSS Thalweg Section", coord_sys=coord_sys) 
    # appoint the section locations 
    locstream["ESMF:Lon"] = lon_thlweg 
    locstream["ESMF:Lat"] = lat_thlweg
    
    gridfile = 'utss_shyfem_esmf_meshinfo.nc'                                 
    srcgrid = ESMF.Mesh(filename=gridfile,filetype=ESMF.FileFormat.ESMFMESH)  
    srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
    
    # salinity
    data = xr.open_dataset(fname, chunks={'node':20000})['salinity']
    ds = data.mean('time')
    secfield = np.zeros((lon_thlweg.shape[0],np.copy(data.level.shape[0])))
    secfield2 = np.zeros((lon_thlweg.shape[0],np.copy(data.level.shape[0])))
    # vertical level kind
    for kind in range(0,93):
        srcfield.data[:] = np.copy(ds[:,kind])
        srcfield.data[srcfield.data==0] = np.NaN
    
        # create a field on the locstream                
        dstfield = ESMF.Field(locstream, name='dstfield')
        dstfield.data[:] = 0.0
    
        # create an object to regrid data from the source to the destination field 
        dst_mask_values=None
        if domask:
                dst_mask_values=np.array([0])
    
        regrid = ESMF.Regrid(srcfield, dstfield,
            # regrid_method=ESMF.RegridMethod.NEAREST_STOD,
            regrid_method=ESMF.RegridMethod.BILINEAR,
            # regrid_method=ESMF.RegridMethod.PATCH,
            unmapped_action=ESMF.UnmappedAction.IGNORE,dst_mask_values=dst_mask_values)
    
        # do the regridding from source to destination field
        dstfield = regrid(srcfield, dstfield)
        secfield[:,kind] = dstfield.data
    
    
    # plt.figure()
    # plt.plot(secfield[0,:],-data.level,'b');
    # plt.plot(tmp_prof[saltcol],-tmp_prof[depthcol],'r');
    
    # salt_utss = np.array([secfield[0,:],-data.level.compute()])
    # salt_obs = np.array([tmp_prof[saltcol],-tmp_prof[depthcol]])
    
    # temperature
    data = xr.open_dataset(fname, chunks={'node':20000})['temperature']
    ds = data.mean('time')
    # vertical level kind
    for kind in range(0,93):
        srcfield.data[:] = np.copy(ds[:,kind])
        srcfield.data[srcfield.data==0] = np.NaN
    
        # create a field on the locstream                
        dstfield = ESMF.Field(locstream, name='dstfield')
        dstfield.data[:] = 0.0
    
        # create an object to regrid data from the source to the destination field 
        dst_mask_values=None
        if domask:
                dst_mask_values=np.array([0])
    
        regrid = ESMF.Regrid(srcfield, dstfield,
            # regrid_method=ESMF.RegridMethod.NEAREST_STOD,
            regrid_method=ESMF.RegridMethod.BILINEAR,
            # regrid_method=ESMF.RegridMethod.PATCH,
            unmapped_action=ESMF.UnmappedAction.IGNORE,dst_mask_values=dst_mask_values)
    
        # do the regridding from source to destination field
        dstfield = regrid(srcfield, dstfield)
        secfield2[:,kind] = dstfield.data
    
    
    # temp_utss = np.array([secfield[0,:],-data.level.compute()])
    # temp_obs = np.array([tmp_prof[tempcol],-tmp_prof[depthcol]])
    
    df1 = pd.DataFrame({"Temp_obs": np.array(tmp_prof[tempcol]),"Salt_obs": np.array(tmp_prof[saltcol]),
                        "zr_obs": -np.array(tmp_prof[depthcol])})
    
    df2 = pd.DataFrame({"Temp_uTSS": secfield2[0,:],"Salt_uTSS": secfield[0,:],
                       "zr_uTSS": -data.level.compute()})
    
    fname = 'obs_ctd_netcdf/2017_01_'+stationame+'_obs.csv'
    df1.to_csv(fname)
    fname = 'obs_ctd_netcdf/2017_01_'+stationame+'_uTSS.csv'
    print(fname)
    df2.to_csv(fname)
    # df2.to_csv(r'obs_ctd/016_09_K0_uTSS.csv')




