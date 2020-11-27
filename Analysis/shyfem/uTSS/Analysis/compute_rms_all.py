import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

root_folder = 'obs_ctd_netcdf/'

stationames = ['K0', '45C', 'MD102', 'M20', 'MD104', 'M14A', 'YSA', 'DIPTAR1Y',
               'DIPTAR2Y', 'DIPTAR3Y', 'DIPTAR4Y', 'MD75', 'AR1', 'MD18',
               'MD101',  'MD13A', 'KD1', 'MD67', 'MD10A', 'MD10B']

# stationames = ['K0', '45C', 'MD102', 'M20', 'MD104', 'M14A', 'YSA', 'DIPTAR1Y',
#                'DIPTAR2Y', 'DIPTAR3Y', 'MD75', 'AR1', 'MD18',
#                 'MD13A', 'KD1', 'MD67', 'MD10A']

# prenames = ['2017_08_', '2018_01_', '2018_04_', '2018_08_']
prenames = ['2018_08_']

# stationames = ['M20', 'M23', 'M8', 'M14', 'MY2', 'MBA', 'M1', 'SB', 'B7',
#                'B13', 'B14']
# prename = '2017_01_'


# plt.figure()
for prename in prenames:
    for ind, stationame in enumerate(stationames):
        print(stationame)
        fname1 = root_folder + prename + stationame + '_bias.csv'
        df1 = pd.read_csv(fname1)
        if ind == 0:
            df = df1**2
        else:
            df = df + df1**2



    rms_salt = np.sqrt(df.Salt_bias/len(stationames))
    rms_temp = np.sqrt(df.Temp_bias/len(stationames))
    df3 = pd.DataFrame({"Temp_rms": np.copy(rms_temp[:len(df1.zr_obs)]),
                        "Salt_rms": np.copy(rms_salt[:len(df1.zr_obs)]),
                       "zr_obs": np.copy(df1['zr_obs'])})

    fname = 'obs_ctd_netcdf/' + prename + 'rms.csv'
    # fname = 'obs_ctd_netcdf_new/' + prename + 'rms.csv'
    df3.to_csv(fname)
