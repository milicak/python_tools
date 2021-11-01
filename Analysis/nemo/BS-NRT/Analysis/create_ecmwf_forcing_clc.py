import numpy as np

ls1 = sorted(glob.glob('/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2014/*/*MED*'))

for ind in np.arange(0,len(ls1)):
    df1 = xr.open_dataset(ls1[ind]) 
    rn_lat1d    =     43.0    #  Column latitude
    rn_lon1d    =     37.0   #  Column longitude
    aa = np.where(np.copy(df1.lon) == rn_lon1d)[0]  
    bb = np.where(np.copy(df1.lat) == rn_lat1d)[0]  
    df1 = df1.isel(lon=aa,lat=bb)
    latitude = np.array([43, 43.125, 43.25])
    longitude = np.array([37, 37.125, 37.25])
    #
    clc = np.tile(df1.TCC,(1,3,3))  
    # clc = clc.transpose(2, 0, 1)    
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
            "CLC": (("time", "latitude", "longitude"), clc),
        },
        {"time": (np.copy(df1.time)), "latitude": latitude, "longitude":
         longitude},
    )
    fname = 'newecmwf_clc_y'+ls1[ind][-22:-18]+'m'+ ls1[ind][-18:-16] + 'd'+ ls1[ind][-16:-14]+'.nc'
    print(fname)
    dd.to_netcdf(fname)
    

# ncks --mk_rec_dmn time forcing_ECMWF_y2016.nc -o dnm.nc
