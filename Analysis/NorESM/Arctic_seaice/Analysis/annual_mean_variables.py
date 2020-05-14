import numpy as np

root_folder = '/cluster/work/users/milicak/RUNS/NorESM/'

member = 'member1'
# member = 'member2'
# member = 'member3'

# variables = ['uflxlvl', 'vflxlvl', 'uhflxlvl', 'vhflxlvl']
# variables = ['TREFHT']
variables = ['mmflxd']
# variables = ['uflxlvl', 'vflxlvl', 'uhflxlvl']

for var in variables:
    for year in range(1850,2015):
        print(var,year)
        # fname = root_folder + member + '/' + member + '.cam.h0.' + var + '.' + str(year) + '*.nc'
        fname = root_folder + member + '/' + member + '.micom.hm.' + var + '.' + str(year) + '*.nc'
        list = sorted(glob.glob(fname))
        df = xr.open_mfdataset(fname)
        time = pd.date_range(np.str(year)+'-01-01', freq='M', periods=12)
        df['time']=time
        ds = df.groupby('time.year').mean('time')
        print(ds.year)
        # ds = df.to_dataset(name=var)
        outname = root_folder + member  + '/' + member + '.micom.hy.'+ var + '.' + str(year) + '.nc'
        # outname = root_folder + member  + '/' + member + '.cam.hy.'+ var + '.' + str(year) + '.nc'
        ds.to_netcdf(outname)






