import glob

def ev_g2c(lambdar, phi, lambda0, phi0):
    ''' dadsasd   '''
    rad = np.pi/180.0
    rEarth = 6378206.4E0
    dlambda = rad * (lambdar - lambda0)                               
    dphi    = rad * (phi - phi0)                                     
    xa = rEarth*dlambda*np.cos(rad*phi0)                                     
    ya = rEarth*dphi      

    return xa,ya

root_folder = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT/'

ls1 = sorted(glob.glob(root_folder + '*nos*.nc'))

df = xr.open_dataset(ls1[0])
df = df.where(df!=0)

dz = np.copy(df.level[1:])-np.copy(df.level[0:-1])
dz = np.append(1,dz)
Nx = df.element.shape[0]
dz2D = np.tile(dz,(Nx,1))
dfs = xr.Dataset({                                      
    'dz': xr.DataArray(                        
                data   = dz2D,                         
                dims   = ['node','level'],                     
        coords = {'level': np.copy(df.level),'node':np.copy(df.element)},  
                ),                                      
})

dlon0 = np.mean(df.longitude)
dlat0 = np.mean(df.latitude)
df['element_index'] -= 1
x1tmp = df.longitude[np.int32(df.element_index[:,0])]
y1tmp = df.latitude[np.int32(df.element_index[:,0])]
x1,y1 = ev_g2c(x1tmp, y1tmp, dlon0, dlat0)
x2tmp = df.longitude[np.int32(df.element_index[:,1])]
y2tmp = df.latitude[np.int32(df.element_index[:,1])]
x2,y2 = ev_g2c(x2tmp, y2tmp, dlon0, dlat0)
x3tmp = df.longitude[np.int32(df.element_index[:,2])]
y3tmp = df.latitude[np.int32(df.element_index[:,2])]
x3,y3 = ev_g2c(x3tmp, y3tmp, dlon0, dlat0)
area = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

salt_time = np.zeros(len(ls1))
for ind, fname in enumerate(ls1):
    print(ind)
    df = xr.open_dataset(fname)
    df = df.where(df!=0)
    df['element_index'] -= 1
    s1tmp = df.salinity[0,np.int32(df.element_index[:,0]),:]
    s2tmp = df.salinity[0,np.int32(df.element_index[:,1]),:]
    s3tmp = df.salinity[0,np.int32(df.element_index[:,2]),:]
    # salt at elements
    salt_el = (s1tmp+s2tmp+s3tmp)/3
    tmp = salt_el.fillna(0)
    mask = xr.where(tmp==0,0,1)
    tmp1 = salt_el*dz2D*mask*area
    tmp2 = dz2D*mask*area 
    salt_time[ind] = tmp1.sum()/tmp2.sum()


time = pd.date_range("2020-01-01", freq="D", periods=len(ls1))
dfs = xr.Dataset({                                                    
    'mean_salt': xr.DataArray(                                             
                data   = salt_time,                                         
                dims   = ['time'],                        
                coords = {'time': time },                      
                attrs  = {                                            
                    'units'     : 'psu'                                
                    }                                                 
                ),                                                    
            },                                                        
    )                                                                 
dfs.to_netcdf('mean_salt.nc')
