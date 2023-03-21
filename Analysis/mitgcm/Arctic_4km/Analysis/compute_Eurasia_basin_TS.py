import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
from matplotlib.path import Path
from scipy.io import loadmat
from pandas.tseries.offsets import DateOffset

# lon1,lat1 is for Kara and Barents Sea
# lon2,lat2 is for Greenland Sea
# lon3,lat3 is for Hudson Bay
# lon4,lat4 is for CAA
# lon5,lat5 is for Arctic Ocean Canadian side
# lon6,lat6 is for Labrador Sea/ Baffin Bay
# lon7,lat7 is for Arctic Ocean Eurasian side
# lon8,lat8 is for Bering Sea
# lon9,lat9 is for Chukchi Sea
# lon10,lat10 is for East Siberian Sea

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid = 'Exp02_0';
# expid = 'Exp02_1';
# expid = 'Exp02_2';
# expid = 'Exp02_3';

datadir = root_folder+expid

ls1 = sorted(glob.glob(datadir+'/2DArcticOcean_1*.nc'))
ls2 = sorted(glob.glob(datadir+'/2DArcticOcean_2*.nc'))

gr = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/Exp02_0/grid.nc')
ls1 = ls1+ls2

df = xr.open_mfdataset(ls1,combine='by_coords')
