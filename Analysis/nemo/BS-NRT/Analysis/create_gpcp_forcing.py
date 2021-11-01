import numpy as np

df = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/sbcatm/bs-nrt_precip-gpcp.nc')
df = df.precip[:,100,390]

latitude = np.array([43, 43.125, 43.25])
longitude = np.array([37, 37.125, 37.25])

precip = np.tile(df,(3,3,1))
precip = precip.transpose(2, 0, 1)    


dd = xr.Dataset(
    {
        "precip": (("time_counter", "latitude", "longitude"), precip),
    },
    {"time_counter": df.time_counter, "latitude": latitude, "longitude":
     longitude},
)

dd.to_netcdf('forcing_gpcp.nc')
