import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
from mpl_toolkits.basemap import Basemap                                            
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

df = pd.read_excel('20190220_1033_ODV_TMAMVerisi.xlsx',sheetname='2016data')

timecol = 'yyyy-mm-ddThh:mm:ss.sss'
depthcol = 'Depth [m]'
tempcol = 'Temperature [degC]'
saltcol = 'Salinity [psu]'

stations = df[timecol].unique()

# for each station
timeind = 0
for station in stations:
    tmp_prof = df.loc[df[timecol]==station]
    if tmp_prof[depthcol].max()>25.0:    
        # compute the day of the year
        year = tmp_prof[timecol][timeind][0:4]
        mnth = tmp_prof[timecol][timeind][5:7]
        day = tmp_prof[timecol][timeind][8:10]
        date = pd.to_datetime(year+mnth+day, format='%Y%m%d')
        new_year_day = pd.Timestamp(year=date.year, month=1, day=1)
        day_of_the_year = (date - new_year_day).days + 1
        day_of_the_year
    
        lon_thlweg = np.array([tmp_prof['Longitude'][timeind]])
        lat_thlweg = np.array([tmp_prof['Latitude'][timeind]])
        
        sdate = "%4.4d" % (day_of_the_year)
        fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
        fname
        
        data = xr.open_dataset(fname, chunks={'node':20000})['salinity']
        ds = data.mean('time')
        
        # dnm = np.vstack([data.longitude,data.latitude])
        # dist = np.zeros(data.longitude.shape)
        
        # np.argmin(dist)
        
        
        domask = False
        # create locstream 
        coord_sys = ESMF.CoordSys.SPH_DEG
        locstream = ESMF.LocStream(lon_thlweg.shape[0], name="uTSS Thalweg Section", coord_sys=coord_sys) 
        # appoint the section locations 
        locstream["ESMF:Lon"] = lon_thlweg 
        locstream["ESMF:Lat"] = lat_thlweg
        
        secfield = np.zeros((lon_thlweg.shape[0],np.copy(data.level.shape[0])))
        gridfile = 'utss_shyfem_esmf_meshinfo.nc'                                 
        srcgrid = ESMF.Mesh(filename=gridfile,filetype=ESMF.FileFormat.ESMFMESH)  
        srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
        
        # vertical level kind
        for kind in range(0,93):
            # print(kind)
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
        

        timeind += tmp_prof.shape[0]
        plt.figure()
        plt.plot(secfield[0,:],-data.level);
        plt.plot(tmp_prof[saltcol],-tmp_prof[depthcol]);
    else:
        timeind += tmp_prof.shape[0]


