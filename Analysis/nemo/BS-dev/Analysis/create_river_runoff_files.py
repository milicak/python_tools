import numpy as np
import glob

root_folder = '/data/opa/bs-mod/upstream_bs-nrt/baseline/sbcrnf/daily/'

ls1 = sorted(glob.glob(root_folder + 'riv*.nc'))

for fname in enumerate(ls1):
    df = xr.open_dataset(fname)
    "ncap2 -s rosaline_20=rosaline_red*0.20 riv_daily_y2014.nc" 
    "ncap2 -s rosaline_25=rosaline_red*0.25 riv_daily_y2014.nc" 
