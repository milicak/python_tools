import numpy as np
from datetime import date
# Load marineHeatWaves definition module
import marineHeatWaves as mhw

df = xr.open_dataset('/work/opa/mi19918/Projects/nemo/Marmara/Marmara_shyfem_SST_2016_2021.nc')

# time for L4 satellite dataset
t = np.arange(date(2016,4,2).toordinal(),date(2021,12,31).toordinal()+1)
dates = [date.fromordinal(tt.astype(int)) for tt in t]



# 17 years times total points
Ntotal = 3654
count = np.zeros((Ntotal,6))
duration = np.zeros((Ntotal,6))
total_days = np.zeros((Ntotal,6))
for ind in range(0,Ntotal):
    print(ind)
    sst = np.copy(df.shyfem_sst[:,ind])
    # mhws, clim = mhw.detect(t, sst)
    mhws, clim = mhw.detect(t, sst,climatologyPeriod=[2020,2021])
    mhwBlock = mhw.blockAverage(t, mhws)
    count[ind] = mhwBlock['count']
    duration[ind] = mhwBlock['duration']
    total_days[ind] = mhwBlock['total_days']



duration[np.isnan(duration)]=0


