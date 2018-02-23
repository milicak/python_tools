from xmitgcm import open_mdsdataset                                             
data_dir = '/export/grunchfs/unibjerknes/milicak/bckup/mitgcm/ice_leads/'   
expid = 'Exp01.3'                                                               
fname = data_dir + expid                                                        
ds = open_mdsdataset(fname,prefix=['S', 'Eta', 'U', 'PTRACER01', 'W', 'V'],ignore_unknown_vars=True)
outname = fname + '/data.nc'
ds.to_netcdf(outname)
