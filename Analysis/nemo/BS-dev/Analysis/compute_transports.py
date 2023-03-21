import numpy as np
import pandas as pd
import os

root_folder = '/work/opa/mi19918/Projects/nemo/BS/'

# expid = 'BS-NRT_MI_bdy7_geom3_rlxbiharm02'
# expid = 'BS-NRT_MI_nemo4_2_02'
expid = 'BS-NRT_MI_bdy7_geom3_rlxbiharm06'

ls1 = glob.glob(root_folder + expid + '/output/volume_transport*')
ls1.sort(key=os.path.getmtime)


df = pd.concat((pd.read_csv(f, header=None, delim_whitespace=True) for f in ls1), ignore_index=True)

# date, time-step, section_number,
# section_name, section_slope_coef, class_number,
# class_name, class_bound_1, class_bound2,
# transport_positive, transport_negative,
# transport_total

# 0        south_sec
# 1        south_sec
# 2       middle_sec
# 3       middle_sec
# 4        north_sec
# 5        north_sec
# 6       after_sill
# 7       after_sill
# 8       west_bound
# 9       west_bound
# 10     west_bound2
# 11     west_bound2
# 12     south_bound
# 13     south_bound
# 14    south_bound2
# 15    south_bound2
# 16      east_bound
# 17      east_bound

# df_new = df.iloc[:, [0,1,9,10,11]]


# north section positive
dsp = df.iloc[5::18,[9]]
# north section negative
dsn = df.iloc[5::18,[10]]
# north section total
dst = df.iloc[5::18,[11]]

# time setup
# time = pd.date_range(str(df.iloc[0,0]), freq="H", periods=df.shape[0])
time = pd.date_range(str(df.iloc[0,0]), freq="8H", periods=dsp.shape[0])

# create dataset
ds = xr.Dataset({
    'Vol_pos': xr.DataArray(
        data   = np.copy(dsp)[:,0],   # enter data here
                dims   = ['time'],
                coords = {'time': time},
                attrs  = {
                    'units'     : 'm3/s'
                    }
                ),
    'Vol_neg': xr.DataArray(
        data   = np.copy(dsn)[:,0],   # enter data here
                dims   = ['time'],
                coords = {'time': time},
                attrs  = {
                    'units'     : 'm3/s'
                    }
                ),
    'Vol_total': xr.DataArray(
        data   = np.copy(dst)[:,0],   # enter data here
                dims   = ['time'],
                coords = {'time': time},
                attrs  = {
                    'units'     : 'm3/s'
                    }
                )
            },
    )

# define data with variable attributes
data_vars = {'vol_positive':(['time'], np.copy(dsp),
                         {'units': 'm3/s'})}

# define coordinates
coords = {'time': (['time'], time)}

# create dataset
ds = xr.Dataset(data_vars=data_vars, 
                coords=coords)
                
