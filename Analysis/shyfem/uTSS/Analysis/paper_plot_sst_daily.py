import numpy as np
from matplotlib import path 

def ev_g2c(lambdar, phi, lambda0, phi0):
    ''' dadsasd   '''
    rad = np.pi/180.0
    rEarth = 6378206.4E0
    dlambda = rad * (lambdar - lambda0)                               
    dphi    = rad * (phi - phi0)                                     
    xa = rEarth*dlambda*np.cos(rad*phi0)                                     
    ya = rEarth*dphi      

    return xa,ya

shapef = {}
def compute_shapefnc(grd):
    dlon0 = 26.68237 
    dlat0 = 40.54126
    x1tmp = grd.longitude[grd.element_index[:,0]]
    y1tmp = grd.latitude[grd.element_index[:,0]]
    x1,y1 = ev_g2c(x1tmp, y1tmp, dlon0, dlat0)
    x2tmp = grd.longitude[grd.element_index[:,1]]
    y2tmp = grd.latitude[grd.element_index[:,1]]
    x2,y2 = ev_g2c(x2tmp, y2tmp, dlon0, dlat0)
    x3tmp = grd.longitude[grd.element_index[:,2]]
    y3tmp = grd.latitude[grd.element_index[:,2]]
    x3,y3 = ev_g2c(x3tmp, y3tmp, dlon0, dlat0)
    area = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    aji = 1.0/area
    b1 = (y2-y3)*aji
    c1 = (x3-x2)*aji
    b2 = (y3-y1)*aji
    c2 = (x1-x3)*aji
    b3 = (y1-y2)*aji
    c3 = (x2-x1)*aji
    dnm = [b1,b2,b3]
    shapef.setdefault('xderiv',[]).append(dnm)
    dnm = [c1,c2,c3]
    shapef.setdefault('yderiv',[]).append(dnm)
    return shapef,area

def crmask(lonlat_polygon, lonlat_model, lonmodel):        
    p = path.Path(lonlat_polygon)                          
    mask = p.contains_points(lonlat_model)                 
    mask.shape = lonmodel.shape                            
    mask = np.multiply(mask, 1)                            
    return mask                                            

df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/daily_clim/uTSS_TS_daily_clim.nc')  
grd = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/uTSS_lobc_chunk_1460.nos.nc')
grd['element_index'] -= 1
shapef,area = compute_shapefnc(grd)

lon = np.copy(df.longitude[0,:])                                    
lat = np.copy(df.latitude[0,:])                                     
lons_lats_model = np.column_stack((lon.flatten(),lat.flatten()))

# marmara sea region
lon1 = np.array([26.83, 27.09, 27.61, 30.06, 30.0, 26.68, 26.53, 26.83]) 
lat1 = np.array([40.55, 40.78, 41.12, 41.06, 40.24, 40.2, 40.41, 40.55])  
mararea = np.column_stack((lon1,lat1))
maskMAR = crmask(mararea,lons_lats_model,lon)                     

sst = df.temperature[:,:,0] 
T1 = sst[:,grd.element_index[:,0]]  
T2 = sst[:,grd.element_index[:,1]]  
T3 = sst[:,grd.element_index[:,2]]  
sst_element = (T1+T2+T3)/3.0 
