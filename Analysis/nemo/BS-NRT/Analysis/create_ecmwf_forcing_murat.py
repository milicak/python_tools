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
ls1 = sorted(glob.glob('/work/opa/mg01720/force/accumulated/IFS_0125/ecmwf_y2014*')) 
ls2 = sorted(glob.glob('/work/opa/mg01720/force/accumulated/IFS_0125/ecmwf_an_y2014*')) 
ls3 = sorted(glob.glob('/work/opa/mg01720/force/accumulated/IFS_0125/ecmwf_spec_y2014*')) 

for ind in np.arange(0,len(ls3)):
    df1 = xr.open_dataset(ls1[ind]) 
    df2 = xr.open_dataset(ls2[ind]) 
    df3 = xr.open_dataset(ls3[ind]) 
    # rn_lat1d    =     43    !  Column latitude
    # rn_lon1d    =     37    !  Column longitude
    df1 = df1.sel(x=112,y=16)
    df2 = df2.sel(x=112,y=16)
    df3 = df3.sel(x=112,y=16)
    latitude = np.array([43, 43.125, 43.25])
    longitude = np.array([37, 37.125, 37.25])
    #
    u10 = np.tile(df1.U10M,(3,3,1))
    u10 = u10.transpose(2, 0, 1)    
    #
    v10 = np.tile(df1.V10M,(3,3,1))
    v10 = v10.transpose(2, 0, 1)    
    #
    t2 = np.tile(df1.T2M,(3,3,1))
    t2 = t2.transpose(2, 0, 1)    
    #
    snw = np.tile(df1.SF,(3,3,1))
    snw = snw.transpose(2, 0, 1)    
    #
    slp = np.tile(df1.SP,(3,3,1))
    slp = slp.transpose(2, 0, 1)    
    #
    precip = np.tile(df2.TP,(3,3,1))
    precip = precip.transpose(2, 0, 1)    
    #
    qsr = np.tile(df2.SSRD,(3,3,1))
    qsr = qsr.transpose(2, 0, 1)    
    #
    qlw = np.tile(df2.STRD,(3,3,1))
    qlw = qlw.transpose(2, 0, 1)    
    #
    rh = np.tile(df3.Q2M,(3,3,1))
    rh = rh.transpose(2, 0, 1)    
    # dd = xr.Dataset(
    #     {
    #         "clc": (("time_counter", "latitude", "longitude"), clc),
    #         "msl": (("time_counter", "latitude", "longitude"), msl),
    #         "precip": (("time_counter", "latitude", "longitude"), precip),
    #         "rh": (("time_counter", "latitude", "longitude"), rh),
    #         "t2": (("time_counter", "latitude", "longitude"), t2),
    #         "u10": (("time_counter", "latitude", "longitude"), u10),
    #         "v10": (("time_counter", "latitude", "longitude"), v10),
    #     },
    #     {"time_counter": df.time_counter, "latitude": latitude, "longitude":
    #      longitude},
    # )
    #
    dd = xr.Dataset(
        {
            "U10m": (("time", "latitude", "longitude"), u10),
            "V10m": (("time", "latitude", "longitude"), v10),
            "SSRD": (("time", "latitude", "longitude"), qsr),
            "STRD": (("time", "latitude", "longitude"), qlw),
            "T2M": (("time", "latitude", "longitude"), t2),
            "Q2M": (("time", "latitude", "longitude"), rh),
            "TP": (("time", "latitude", "longitude"), precip),
            "SF": (("time", "latitude", "longitude"), snw),
            "SP": (("time", "latitude", "longitude"), slp),
        },
        {"time": (np.copy(df1.time)), "latitude": latitude, "longitude":
         longitude},
    )
    fname = 'new'+ls1[ind][-20:] 
    print(fname)
    dd.to_netcdf(fname)
    

# ncks --mk_rec_dmn time forcing_ECMWF_y2016.nc -o dnm.nc
