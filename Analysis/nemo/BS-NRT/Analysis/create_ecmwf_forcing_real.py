import numpy as np

#   sn_wndi     = 'ecmwf'                 ,    6              , 'U10M'    ,   .true.   , .false. , 'daily'  , 'weights_bicub.nc'            , ''       , ''
#   sn_wndj     = 'ecmwf'                 ,    6              , 'V10M'    ,   .true.   , .false. , 'daily'  , 'weights_bicub.nc'            , ''       , ''
#   sn_qsr      = 'ecmwf_an'              ,    6              , 'SSRD'    ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''
#   sn_qlw      = 'ecmwf_an'              ,    6              , 'STRD'    ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''
#   sn_tair     = 'ecmwf'                 ,    6              , 'T2M'     ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''
#   sn_humi     = 'ecmwf_spec'            ,    6              , 'Q2M'     ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''
#   sn_prec     = 'ecmwf_an'              ,    6              , 'TP'      ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''
#   sn_snow     = 'ecmwf'                 ,    6              , 'SF'      ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''
#   sn_slp      = 'ecmwf'                 ,    6              , 'SP'      ,   .true.   , .false. , 'daily'  , 'weights_bilin.nc'            , ''       , ''

gr = xr.open_dataset('/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2014/12/20141201-ECMWF---AM0125-BLK-b20141202_an-fv01.00.nc')

ls1 = sorted(glob.glob('/data/opa/bs-mod/upstream_bs-nrt/sbcatm2/*y2016*'))
ls2 = sorted(glob.glob('/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2016/*/*MED*fv06*'))

for ind in np.arange(0,len(ls1)):
    df1 = xr.open_dataset(ls1[ind]) 
    df2 = xr.open_dataset(ls2[ind]) 
    # rn_lat1d    =     43    !  Column latitude
    # rn_lon1d    =     37    !  Column longitude
    rn_lat1d    =     43.0    #  Column latitude
    rn_lon1d    =     37.0   #  Column longitude
    aa = np.where(np.copy(gr.lon) == rn_lon1d)[0]  
    bb = np.where(np.copy(gr.lat) == rn_lat1d)[0]  
    df1 = df1.sel(x=aa,y=bb)
    aa = np.where(np.copy(df2.lon) == rn_lon1d)[0]  
    bb = np.where(np.copy(df2.lat) == rn_lat1d)[0]  
    df2 = df2.isel(lon=aa,lat=bb)
    latitude = np.array([43, 43.125, 43.25])
    longitude = np.array([37, 37.125, 37.25])
    #
    u10 = np.tile(df1.u10,(1,3,3))
    #
    v10 = np.tile(df1.v10,(1,3,3))
    #
    t2 = np.tile(df1.t2,(1,3,3))
    #
    msl = np.tile(df1.msl,(1,3,3))*1e-2
    #
    rh = np.tile(df1.rh,(1,3,3))
    #
    q2m = np.tile(df2.D2M,(1,3,3))
    #
    ssrd = np.tile(df1.ssrd,(1,3,3))
    #
    strd = np.tile(df1.strd,(1,3,3))
    #
    sf = np.tile(df1.sf,(1,3,3))
    #
    tp = np.tile(df1.tp,(1,3,3))
    #
    clc = np.tile(df2.TCC,(1,3,3))
    #
    dd = xr.Dataset(
        {
            "u10": (("time", "latitude", "longitude"), u10),
            "v10": (("time", "latitude", "longitude"), v10),
            "t2": (("time", "latitude", "longitude"), t2),
            "rh": (("time", "latitude", "longitude"), rh),
            "q2m": (("time", "latitude", "longitude"), q2m),
            "msl": (("time", "latitude", "longitude"), msl),
            "ssrd": (("time", "latitude", "longitude"), ssrd),
            "strd": (("time", "latitude", "longitude"), strd),
            "sf": (("time", "latitude", "longitude"), sf),
            "tp": (("time", "latitude", "longitude"), tp),
            "clc": (("time", "latitude", "longitude"), clc),
        },
        {"time": (np.copy(df1.time_counter)), "latitude": latitude, "longitude":
         longitude},
    )
    fname = 'new'+ls1[ind][-20:] 
    print(fname)
    dd.to_netcdf(fname)
    

# ncks --mk_rec_dmn time forcing_ECMWF_y2016.nc -o dnm.nc
# time = pd.date_range("2016-01-01", freq="6H", periods=366 * 4) 
