import numpy as np

df = xr.open_dataset('forcing_ECMWF_y2016_org.nc')

latitude = np.array([43, 43.125, 43.25])
longitude = np.array([37, 37.125, 37.25])

clc = np.tile(df.clc,(3,3,1))
clc = clc.transpose(2, 0, 1)    

msl = np.tile(df.msl,(3,3,1))
msl = msl.transpose(2, 0, 1)    

precip = np.tile(df.precip,(3,3,1))
precip = precip.transpose(2, 0, 1)    

rh = np.tile(df.rh,(3,3,1))
rh = rh.transpose(2, 0, 1)    

t2 = np.tile(df.t2,(3,3,1))
t2 = t2.transpose(2, 0, 1)    

u10 = np.tile(df.u10,(3,3,1))
u10 = u10.transpose(2, 0, 1)    

v10 = np.tile(df.v10,(3,3,1))
v10 = v10.transpose(2, 0, 1)    

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
        "clc": (("time", "latitude", "longitude"), clc),
        "msl": (("time", "latitude", "longitude"), msl),
        "precip": (("time", "latitude", "longitude"), precip),
        "rh": (("time", "latitude", "longitude"), rh),
        "t2": (("time", "latitude", "longitude"), t2),
        "u10": (("time", "latitude", "longitude"), u10),
        "v10": (("time", "latitude", "longitude"), v10),
    },
    {"time": np.float(np.copy(df.time_counter)), "latitude": latitude, "longitude":
     longitude},
)

dd.to_netcdf('forcing_ECMWF_y2016.nc')

# ncks --mk_rec_dmn time forcing_ECMWF_y2016.nc -o dnm.nc
