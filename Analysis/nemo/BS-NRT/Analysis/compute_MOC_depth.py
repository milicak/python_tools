import numpy as np


root_folder = '/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/reanalysis/day/'

# gr = xr.open_dataset('/data/opa/bs-mod/experiments/bs-simu_6.6_zeus/36_1d_20140903_20140909_grid_T.nc')
gr = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/mesh_mask_Aug19_31-5-25.nc')
gr = gr.rename({'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'})
# year = 1993
# files = root_folder + np.str(year) + '/*/*RFVL*BSe3r1*'      
files = root_folder + '*/*/*RFVL*BSe3r1*'      
ls1 = sorted(glob.glob(files)) 
# ind = 0
for ind in range(0,len(ls1)):
    df = xr.open_dataset(ls1[ind])  
    # zonal transport                                           
    voltrU = df.uo*gr.e3u_0*gr.e2u 
    # voltrU = df.uo*gr.e3u_0[0,:,:]*gr.e2u[0,:,:]
    Try = voltrU.sum(('lat'))                                     
    # Trymean = Try.mean('time')                                  
    # reverse z coordinate for different cumsum                 
    Tryreverse_z = Try.reindex(depth=Try.depth[::-1])       
    # Trymeansum = Trymean.cumsum('depth')                                                
    Tryreverse_zsum = Tryreverse_z.cumsum('depth')                              
    Tryreverse_zsum = Tryreverse_zsum.reindex(depth=Tryreverse_zsum.depth[::-1])
    Tryreverse_zsum = Tryreverse_zsum.to_dataset(name='MOC_zonal')
    print(ls1[ind][79:87])
    outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/moc_depth_zonal_BSe3r1_' + ls1[ind][79:87] + '.nc'
    Tryreverse_zsum.to_netcdf(outname)
