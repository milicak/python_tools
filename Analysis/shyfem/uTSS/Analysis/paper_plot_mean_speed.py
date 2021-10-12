import numpy as np
from scipy.spatial import cKDTree
from mpl_toolkits.basemap import Basemap
import xarray as xr
import cmocean

def speedToLW(speed):
    ''' 
    Function to convert windspeed into a sensible linewidth
    This will need to change depending on your data
    '''
    return 0.5 + 3.5 * speed 

def lon_lat_to_cartesian(lon, lat, R = 6371000):                                
    """                                                                         
    Calculates lon, lat coordinates of a point on a sphere with                 
    radius R. Taken from http://earthpy.org/interpolation_between_grids_with_ckdtree.html
                                                                                
    Parameters                                                                  
    ----------                                                                  
    lon : 1d array                                                              
        longitudes                                                              
    lat : 1d array                                                              
        latitudes                                                               
    R   : float                                                                 
        radius of the sphere                                                    
                                                                                
    Returns                                                                     
    -------                                                                     
    x,y,z : 1d arrays                                                           
        cartesian coordinates                                                   
    """                                                                         
    lon_r = np.radians(lon)                                                     
    lat_r = np.radians(lat)                                                     
                                                                                
    x =  R * np.cos(lat_r) * np.cos(lon_r)                                      
    y =  R * np.cos(lat_r) * np.sin(lon_r)                                      
    z =  R * np.sin(lat_r)                                                      
    return x,y,z                                                                


def create_indexes_and_distances(lon_gr, lat_gr, lons, lats, kint=1, n_jobs=2, ):            
    '''                                                                         
    Creates KDTree object and query it for indexes of points in FESOM mesh that are close to the
    points of the target grid. Also return distances of the original points to target points.
                                                                                
                                                                                
    Parameters                                                                  
    ----------                                                                  
    lon_gr : longititude shyfem                                                 
        pyfesom mesh representation                                             
    lat_gr :  latitude shyfem
    lons/lats : array                                                           
        2d arrays with target grid values.                                      
    k : int                                                                     
        k-th nearest neighbors to return.                                       
    n_jobs : int, optional                                                      
        Number of jobs to schedule for parallel processing. If -1 is given      
        all processors are used. Default: 1.                                    
                                                                                
    Returns                                                                     
    -------                                                                     
    distances : array of floats                                                 
        The distances to the nearest neighbors.                                 
    inds : ndarray of ints                                                      
        The locations of the neighbors in data.                                 
                                                                                
    '''                                                                         
    xs, ys, zs = lon_lat_to_cartesian(lon_gr, lat_gr)                         
    xt, yt, zt = lon_lat_to_cartesian(lons.flatten(), lats.flatten())           
                                                                                
    tree = cKDTree(list(zip(xs, ys, zs)))                                       
    distances, inds = tree.query(list(zip(xt, yt, zt)), k = kint, n_jobs=n_jobs)   
                                                                                
    return distances, inds                                                      

def shyfem2regular(data, lon_gr, lat_gr, lons, lats, distances=None,                      
                  inds=None, how='nn', k=1, radius_of_influence=100000, n_jobs = 2 ):
    '''                                                                         
    Interpolates data from Shyfem mesh to target (usually regular) mesh.         
                                                                                
    Parameters                                                                  
    ----------                                                                  
    data : array                                                                
        1d array that represents FESOM data at one level.                       
    mesh : fesom_mesh object                                                    
        pyfesom mesh representation                                             
    lons/lats : array                                                           
        2d arrays with target grid values.                                      
    distances : array of floats, optional                                       
        The distances to the nearest neighbors.                                 
    inds : ndarray of ints, optional                                            
        The locations of the neighbors in data.                                 
    how : str                                                                   
       Interpolation method. Options are 'nn' (nearest neighbor) and 'idist' (inverce distance)
    k : int                                                                     
        k-th nearest neighbors to use. Only used when how==idist                
    radius_of_influence : int                                                   
        Cut off distance in meters.                                             
    n_jobs : int, optional                                                      
        Number of jobs to schedule for parallel processing. If -1 is given      
        all processors are used. Default: 1.                                    
                                                                                
    Returns                                                                     
    -------                                                                     
    data_interpolated : 2d array                                                
        array with data interpolated to the target grid.                        
                                                                                
    '''                                                                         
     #print distances                                                            
    if (distances is None) or (inds is None):                                   
                                                                                
        if how=='nn':                                                           
            distances, inds = create_indexes_and_distances(lon_gr, lat_gr,
                                                           lons, lats, kint=1, n_jobs=n_jobs)  
        elif how=='idist':                                                      
            distances, inds = create_indexes_and_distances(lon_gr, lat_gr,
                                                           lons, lats, kint=k, n_jobs=n_jobs)  
                                                                                
    if distances.ndim == 1:                                                     
        #distances_ma = np.ma.masked_greater(distances, radius_of_influence)    
        data_interpolated = data[inds]                                          
                                                                                
        data_interpolated[distances>=radius_of_influence] = np.nan              
                                                                                
        data_interpolated = data_interpolated.reshape(lons.shape)               
        data_interpolated = np.ma.masked_invalid(data_interpolated)             
    else:                                                                       
        distances_ma = np.ma.masked_greater(distances, radius_of_influence)     
                                                                                
        w = 1.0 / distances_ma**2                                               
        data_interpolated = np.ma.sum(w * data[inds], axis=1) / np.ma.sum(w, axis=1)
        data_interpolated.shape = lons.shape                                    
        data_interpolated = np.ma.masked_invalid(data_interpolated)             
                                                                                
    return data_interpolated                                                    


root_folder = '/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/'

fname = root_folder + 'uTSS_UV_monthly_clim.nc'

df = xr.open_dataset(fname)

du = df.u_velocity.mean('month')
dv = df.v_velocity.mean('month')

sp = (du[:,0]**2)+(dv[:,0]**2) 

m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

uvel = np.copy(du[:,0])
vvel = np.copy(dv[:,0])
# regular grid 
# x = np.linspace(22.5, 32.0, 500)
# y = np.linspace(38.5, 43.0, 400)
# x=np.arange(25.0,30.0,0.01)
# y=np.arange(39.5,42.0,0.01)
# xi, yi = np.meshgrid(x,y)

lons, lats, x, y = m.makegrid(100, 100, returnxy=True) 

u_int = shyfem2regular(uvel, np.copy(df.longitude[0,:]),
                       np.copy(df.latitude[0,:]), lons, lats,
                       radius_of_influence=50000)
v_int = shyfem2regular(vvel, np.copy(df.longitude[0,:]),
                       np.copy(df.latitude[0,:]), lons, lats,
                       radius_of_influence=50000)

plt.figure()
m.drawcoastlines(linewidth=0.2);
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(df.longitude[0,:]),np.copy(df.latitude[0,:]))

im1 = plt.tripcolor(longitude,latitude,df.element_index[-1,:,:]-1,np.sqrt(sp)
             ,vmin=0,vmax=0.6,shading='gouraud',cmap=cmocean.cm.speed);

speed = np.copy(np.sqrt(u_int**2+v_int**2))
speed = speed/speed.max()
m.streamplot(x, y, u_int, v_int, color='k', density=[3,2], linewidth=speedToLW(speed))
m.fillcontinents(color='grey',zorder=2);

cb = m.colorbar(im1,"right", size="5%", pad="4%", extend='max')
cb.set_label('[m/s]',rotation=0,y=1.07,labelpad=-30) 
plt.savefig('paperfigs/uTSS_mean_speed_streamplot.png', bbox_inches='tight',format='png',dpi=300)

