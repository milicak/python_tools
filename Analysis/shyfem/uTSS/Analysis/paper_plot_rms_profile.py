import numpy as np
import xarray as xr

root_folder1 = 'obs_ctd_netcdf/'
root_folder2 = 'obs_ctd_netcdf_athena/'


prenames = ['2017_08_', '2018_01_', '2018_04_', '2018_08_']


fnames = ['obs_ctd_netcdf_athena/2017_08_',
          'obs_ctd_netcdf/2018_01_',
          'obs_ctd_netcdf_athena/2018_04_',
          'obs_ctd_netcdf_athena/2018_08_']

for fname in fnames:
    fname = fname + 'rms.csv'
    df = pd.read_csv(fname)
    plt.plot(df.Salt_rms*0.8,df.zr_obs);


plt.legend(['2017-08','2018-01','2018-04','2018-08'])
plt.xlim(0,5);
plt.savefig('paperfigs/salinity_rms_all.png', bbox_inches='tight',format='png',dpi=300)


for fname in fnames:
    fname = fname + 'rms.csv'
    df = pd.read_csv(fname)
    plt.plot(df.Temp_rms*0.8,df.zr_obs);


plt.legend(['2017-08','2018-01','2018-04','2018-08'])
plt.xlim(0,5);
plt.savefig('paperfigs/temperature_rms_all.png', bbox_inches='tight',format='png',dpi=300)
