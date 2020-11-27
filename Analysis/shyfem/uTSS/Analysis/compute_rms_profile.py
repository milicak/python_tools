import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

root_folder = 'obs_ctd_netcdf/'

stationames = ['K0', '45C', 'MD102', 'M20', 'MD104', 'M14A', 'YSA', 'DIPTAR1Y',
               'DIPTAR2Y', 'DIPTAR3Y', 'DIPTAR4Y', 'MD75', 'AR1', 'MD18',
               'MD101',  'MD13A', 'KD1', 'MD67', 'MD10A', 'MD10B']

prename = '2018_08_'

# stationames = ['M20', 'M23', 'M8', 'M14', 'MY2', 'MBA', 'M1', 'SB', 'B7',
#                'B13', 'B14']
# prename = '2017_01_'

for stationame in stationames:
    print(stationame)
    fname1 = root_folder + prename + stationame + '_obs.csv'
    fname2 = root_folder + prename + stationame + '_uTSS.csv'
    df1 = pd.read_csv(fname1)
    df2 = pd.read_csv(fname2)
    ind_last = np.sum(~np.isnan(df2['Salt_uTSS']))
    df2['Salt_uTSS'][ind_last:] = df2['Salt_uTSS'][ind_last-1]
    df2['Temp_uTSS'][ind_last:] = df2['Temp_uTSS'][ind_last-1]
    # f_nearest = interp1d(df2['zr_uTSS'], df2['Salt_uTSS'], kind='nearest')
    f_linear  = interp1d(df2['zr_uTSS'], df2['Salt_uTSS'])
    salt_int = f_linear(df1['zr_obs'])
    f_linear  = interp1d(df2['zr_uTSS'], df2['Temp_uTSS'])
    temp_int = f_linear(df1['zr_obs'])
    # compute biases
    salt_bias = salt_int - df1['Salt_obs']
    temp_bias = temp_int - df1['Temp_obs']
    df3 = pd.DataFrame({"Temp_bias": np.copy(temp_bias),"Salt_bias":
                        np.copy(salt_bias),
                       "zr_obs": np.copy(df1['zr_obs'])})

    fname = 'obs_ctd_netcdf/' + prename + stationame + '_bias.csv'
    # fname = 'obs_ctd_netcdf_new/' + prename + stationame + '_bias.csv'
    df3.to_csv(fname)





