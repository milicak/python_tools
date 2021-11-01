import numpy as np

# from 1993 - 2019
root_folder = '/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/reanalysis/day/'
# for 2020
root_folder = '/data/products/BSFS/bs-rea_v3.1/cmems/reanalysis_daily_mean/'


# gr = xr.open_dataset('/data/opa/bs-mod/experiments/bs-simu_6.6_zeus/36_1d_20140903_20140909_grid_T.nc')
gr = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/mesh_mask_Aug19_31-5-25.nc')
gr = gr.rename({'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'})
# year = 1993
# files = root_folder + np.str(year) + '/*/*RFVL*BSe3r1*'      
# from 1993 - 2019
files = root_folder + '*/*/*RFVL*BSe3r1*'      
# for 2020
files = root_folder + np.str(2020) +'/*RFVL*BSe3r1*'      
ls1 = sorted(glob.glob(files)) 
# ind = 0
for ind in range(0,len(ls1)):
    df = xr.open_dataset(ls1[ind])  
    # meridional transport                                           
    voltrV = df.vo*gr.e3v_0*gr.e1v 
    Trx = voltrV.sum(('lon'))                                     
    # reverse z coordinate for different cumsum                 
    Trxreverse_z = Trx.reindex(depth=Trx.depth[::-1])       
    Trxreverse_zsum = Trxreverse_z.cumsum('depth')                              
    Trxreverse_zsum = Trxreverse_zsum.reindex(depth=Trxreverse_zsum.depth[::-1])
    Trxreverse_zsum = -Trxreverse_zsum
    Trxreverse_zsum = Trxreverse_zsum.to_dataset(name='MOC_meridional')
    # from 1993 - 2019 
    # print(ls1[ind][79:87])
    # outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/moc_depth_meridional_BSe3r1_' + ls1[ind][79:87] + '.nc'
    # for 2020
    print(ls1[ind][65:73])
    outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/moc_depth_meridional_BSe3r1_' + ls1[ind][65:73] + '.nc'
    Trxreverse_zsum.to_netcdf(outname)
