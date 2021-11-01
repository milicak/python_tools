import numpy as np

def coriolis(lat):
    """Compute the Coriolis parameter for the given latitude:
    ``f = 2*omega*sin(lat)``, where omega is the angular velocity 
    of the Earth.
    
    Parameters
    ----------
    lat : array
      Latitude [degrees].
    """
    omega   = 7.2921159e-05  # angular velocity of the Earth [rad/s]
    return 2*omega*np.sin(lat/360.*2*np.pi)


gr1 = xr.open_dataset('/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_2.0/model_old/domain_cfg_flip_cor.nc')

# newlat = np.linspace(27.25,42,591)   
newlat = np.linspace(27.25-18,60,591) 
dx = np.copy(newlat[1])-np.copy(newlat[0]) 
newlat_f = np.linspace(27.25-18+0.5*dx,60+0.5*dx,591) 

newlat = np.tile(newlat,(261,1)) 
newfcor_t = coriolis(newlat)
newfcor_t = np.fliplr(newfcor_t)

newlat_f = np.tile(newlat_f,(261,1)) 
newfcor_f = coriolis(newlat_f)
newfcor_f = np.fliplr(newfcor_f)

gr1['ff_f'][0,:,:] = newfcor_f                                                                                                                                        
gr1['ff_t'][0,:,:] = newfcor_t   

gr1.to_netcdf('/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_2.0/model_old/domain_cfg_flip_corv2.nc')
