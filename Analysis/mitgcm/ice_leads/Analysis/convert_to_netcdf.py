from xmitgcm import open_mdsdataset
import os

data_dir = '/export/grunchfs/unibjerknes/milicak/bckup/mitgcm/ice_leads/'
expid = 'Exp01.12'
fname = data_dir + expid
ds = open_mdsdataset(fname,prefix=['S', 'Eta', 'U', 'HEFF', 'HSNOW', 'PTRACER01', 'W', 'V'],ignore_unknown_vars=True)
outname = fname + '/data.nc'
ds.to_netcdf(outname)

ds = open_mdsdataset(fname,prefix=['Qnet', 'Qsw', 'EmPmR', 'FU', 'FV'],ignore_unknown_vars=True)
outname = fname + '/data2.nc'
ds.to_netcdf(outname)

# rm Qnet*.meta Qnet*.data EmPmR*.meta EmPmR*.data Qsw*.meta Qsw*.data FU*.meta FU*.data FV*.meta FV*.data
# rm S*.meta S*.data Eta*.meta Eta*.data U*.meta U*.data V*.meta V*.data
# rm W*.meta W*.data PTRACER01*.meta PTRACER01*.data
# rm HEFF*.meta HEFF*.data HSNOW*.meta HSNOW*.data

#rmname = fname + '/S*meta'
#cmmnd = 'rm ' + rmname
#os.system(cmmnd)
