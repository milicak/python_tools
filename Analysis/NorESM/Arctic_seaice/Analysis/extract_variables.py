import numpy as np

root_folder = '/tos-project3/NS9560K/noresm/cases/'
out_folder = '/cluster/work/users/milicak/RUNS/NorESM/'

member = 'member1'
# member = 'member2'
# member = 'member3'

# variables = ['fice', 'uflxlvl', 'vflxlvl', 'uhflxlvl', 'vhflxlvl']
# variables = ['TREFHT']
variables = ['mmflxd']

expid = 'NHIST_f19_tn14_20190625'
# expid = 'NHIST_02_f19_tn14_20190801'
# expid = 'NHIST_03_f19_tn14_20190801'
for var in variables:
    for year in range(1850,1950):
        for mnth in range(1,13):
            print(var,year,mnth)
            fname = root_folder + expid + '/ocn/hist/' + expid + '.micom.hm.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            # fname = root_folder + expid + '/atm/hist/' + expid + '.cam.h0.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            df = xr.open_dataset(fname)[var]
            ds = df.to_dataset(name=var)
            outname = out_folder + member  + '/' + member + '.micom.hm.'+ var + '.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            # outname = out_folder + member  + '/' + member + '.cam.h0.'+ var + '.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            ds.to_netcdf(outname)






expid = 'NHIST_f19_tn14_20190710'
# expid = 'NHIST_02_f19_tn14_20190813'
# expid = 'NHIST_03_f19_tn14_20190813'
for var in variables:
    for year in range(1950,2015):
        for mnth in range(1,13):
            print(var,year,mnth)
            fname = root_folder + expid + '/ocn/hist/' + expid + '.micom.hm.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            # fname = root_folder + expid + '/atm/hist/' + expid + '.cam.h0.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            df = xr.open_dataset(fname)[var]
            ds = df.to_dataset(name=var)
            outname = out_folder + member  + '/' + member + '.micom.hm.'+ var + '.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            # outname = out_folder + member  + '/' + member + '.cam.h0.'+ var + '.' + str(year) + '-' + str(mnth).zfill(2) + '.nc'
            ds.to_netcdf(outname)




