from datetime import date, timedelta
from PyCNAL_regridding import *
import glob

def month_range(start, periods=12):
    rng = pd.date_range(pd.Timestamp(start)-pd.offsets.MonthBegin(),
                        periods=periods,
                        freq='MS')
    ret = (rng + pd.offsets.Day(pd.Timestamp(start).day-1)).to_series()
    ret.loc[ret.dt.month > rng.month] -= pd.offsets.MonthEnd(1)
    return pd.DatetimeIndex(ret)

def remove_extreme_val1(ds,var,val1):
    ds[var] = ds[var].where(ds[var]>val1)
    return ds

def remove_extreme_val1_val2(ds,var,val1,val2):
    ds[var] = ds[var].where((ds[var]<val1) & (ds[var]>val2))
    return ds


daystart = 21184
time = month_range('1958-01-15', 12*61)
time = xr.cftime_range(start="1958-01-01", periods=12*61, freq="1MS", calendar="gregorian")
time2 = time.shift(15,'D')

# OMIP2.025d.01_1m_19720101_19721231_grid_T_obc_10.nc
years=np.arange(1958,1959)
# years=np.arange(1958,1982)
# years=np.arange(1982,1982)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
vars=('temp_segment_001','salt_segment_001','u_segment_001','v_segment_001',
      'temp_segment_002','salt_segment_002','u_segment_002','v_segment_002',
      'temp_segment_003','salt_segment_003','u_segment_003','v_segment_003',
      'temp_segment_004','salt_segment_004','u_segment_004','v_segment_004')

for ind, year in enumerate(years):
    print(year)
    ls1 = sorted(glob.glob(mom_dir+'*'+str(year)+'*grid_T_obc*'))
    df = xr.open_mfdataset(ls1, decode_times=False)
    # df['temp_segment_001'] = df.temp_segment_001.where(df.temp_segment_001>-1.8)
    # df['u_segment_001'] = df.u_segment_001.where((df.u_segment_001<2) & (df.u_segment_001>-2))
    for var in vars:
        print(var)
        if var[0] == 't':
            remove_extreme_val1(df,var,-1.8);
            df[var] = df[var].ffill(dim=df[var].dims[1])
            df[var] = df[var].fillna(0)
        elif var[0] == 's':
            remove_extreme_val1(df,var,0);
            df[var] = df[var].ffill(dim=df[var].dims[1])
        elif var[0] == 'u':
            remove_extreme_val1_val2(df,var,1,-1);
            df[var] = df[var].fillna(0)
            # df[var][:,28::,:,:]=0
        elif var[0] == 'v':
            remove_extreme_val1_val2(df,var,1,-1);
            df[var] = df[var].fillna(0)
            # df[var][:,28::,:,:]=0
    
    encoding={'time':{'dtype':'float64'}}
    # df.time.attrs['units'] = 'Seconds since 01/01/1958 00:00:00 UTC'
    # df.time.attrs['calendar']='gregorian'
    df['time'] = time[ind*12:(12*(ind+1))]
    fout = mom_dir + 'MOM6_Arctic_Copernicus_year_' + str(year) + '_obc.nc'
    df.to_netcdf(fout, encoding=encoding, unlimited_dims={'time':True})
    # df['time'] = df.time+365*ind
    # df['time'].attrs['calendar']='gregorian'

# d0 = date(1900, 1, 1)
# days = np.array([15,14,15,15,15,15,15,15,15,15,15,15])
# aa=xr.cftime_range(start="1958-01-14", periods=1, freq="1MS", calendar="gregorian")
# for ind in range(0,len(ls1)):
#     df = xr.open_dataset(ls1[ind])
#     year = int(ls1[ind][-25:-21])
#     month = int(ls1[ind][-5:-3])+1
#     if month==2:
#         days=14
#     else:
#         days=15
#
#     d1 = date(year, month, days)
#     delta = d1 - d0
#     days = delta.days
#     print(days)
#     time = timeobject(days)
#     time.units = 'days since 1900-01-01'
#     time.calendar = 'gregorian'
