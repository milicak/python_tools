import numpy as np


def nsidcmn(fname):
    ns = pd.read_csv(fname)
    # extract 1992 to 2016
    ns = ns[13:38]
    mnmean = ns[' extent'].mean()
    return mnmean

nsidc_ext = np.array([12.11,
    11.923,
    12.011,
    11.415,
    11.841,
    11.668,
    11.757,
    11.691,
    11.508,
    11.6,
    11.363,
    11.397,
    11.24,
    10.907,
    10.773,
    10.474,
    10.978,
    10.932,
    10.711,
    10.483,
    10.406,
    10.897,
    10.79,
    10.566,
    10.163,
    10.393,
    10.355])

nsidc_area = np.array([10.086,
    9.786,
    9.968,
    9.403,
    9.791,
    9.593,
    9.601,
    9.672,
    9.445,
    9.545,
    9.332,
    9.265,
    9.356,
    8.998,
    8.779,
    8.434,
    9.221,
    9.21,
    8.955,
    8.749,
    8.65,
    9.215,
    9.085,
    8.895,
    8.311,
    8.684,
    8.732])

nsidc_year = np.arange(1992,2019)

fname = 'Exp02_0_seaice_extent.nc'
df = xr.open_dataset(fname)
siextann = df.SI_extent.groupby('time.year').mean('time') 

fname = 'Exp02_0_seaice_area.nc'
df2 = xr.open_dataset(fname)
siareaann = df2.SIarea.groupby('time.year').mean('time') 

siextmth = df.SI_extent.groupby('time.month').mean('time')   
siareamth = df2.SIarea.groupby('time.month').mean('time')   

nsimean = np.zeros(12)
for ind in range(0,12):
    fname = 'nsidc/N_'+ str(ind+1).zfill(2) +'_extent_v3.0.csv' 
    nsimean[ind] = nsidcmn(fname)





