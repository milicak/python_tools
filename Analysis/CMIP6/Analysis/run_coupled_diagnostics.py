import xarray as xr
import numpy as np
import glob
from datetime import datetime
from dask.distributed import Client, LocalCluster
#                                                                                                                                                                                                                                                                              
from OMIP_utils import annual_mean, annual_mean_loop, thermosteric_sealvl, thermosteric_sealvl_ref_period, upper_700m_mean

if __name__ == '__main__':
    from dask.distributed import Client, LocalCluster                                                                                                                                                                                                                          
    cluster = LocalCluster(n_workers=6, memory_limit=25e9, local_dir='/cluster/projects/nn2345k/documentation_figures_BLOM/diagnostics_coupled/NorESM2-MM/output/')                                                                                                            
    client = Client(cluster)        

machine='FRAM'
for expID in ['piControl','historical1']:
    print(machine,expID)
    if expID in ['piControl']:
        if machine in ['NIRD']:   
            fpath = '/tos-project3/NS9560K/noresm/cases/' #data on NIRD
            case  = 'N1850frc2_f09_tn14_20191001'
        elif machine in ['FRAM']:
            fpath = '/cluster/NS9560K/noresm/cases/'      #data on FRAM
            case  = 'N1850frc2_f09_tn14_20191012'
    elif expID in ['historical1']:
        if machine in ['NIRD']:
            fpath = '/tos-project3/NS9560K/noresm/cases/' #data on NIRD
            case  = 'NHISTfrc2_f09_tn14_20191001'
        elif  machine in['FRAM']:                          #data on FRAM
            fpath = '/cluster/NS9560K/noresm/cases/'
            case  = 'NHISTfrc2_f09_tn14_20191025'

    # 
    fnames_ocn = sorted(glob.glob(fpath+case+'/ocn/hist/*.micom.hm*.nc'))
    fnames_ice = sorted(glob.glob(fpath+case+'/ice/hist/*cice.h.????-??.nc'))
    # JUST THE FULL PERIOD
    year_range1=slice(0,len(fnames_ocn))
    year_range1_ice=slice(0,len(fnames_ice))
    #
    exec(open('create_coupled_diagnostics.py').read())
