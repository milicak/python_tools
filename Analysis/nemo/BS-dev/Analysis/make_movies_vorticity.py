import numpy as np 
import cartopy.feature                                                                
import cartopy.crs as ccrs                                                            
                                                                                      

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

def rel_vort(dfu, dfv, gr):
    uvel = dfu.uo[:,0,:,:]*gr.e1u                
    vvel = dfv.vo[:,0,:,:]*gr.e2v                
    du = uvel[:,1:,:]-uvel[:,:-1,:]                  
    du = du.to_dataset(name='du')                
    tmp = du.du[:,-1,:]                            
    tmp = tmp.to_dataset(name='du')              
    du = xr.concat([du,tmp],dim='y')             
    dv = vvel[:,:,1:]-vvel[:,:,:-1]                  
    dv = dv.to_dataset(name='dv')                
    tmp = dv.dv[:,:,-1]                            
    tmp = tmp.to_dataset(name='dv')              
    dv = xr.concat([dv,tmp],dim='x')             
    zeta = 1/(gr.e1f*gr.e2f) * (dv.dv - du.du)   
    return zeta

gr = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')
gr = gr.isel(t=0)                            
# compute fcor
fcor = coriolis(gr.nav_lat)

root_folder = '/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_nemo4_2_07/rebuilt/'

ls1t = sorted(glob.glob(root_folder+'*grid_T.nc'))
ls1u = sorted(glob.glob(root_folder+'*grid_U.nc'))
ls1v = sorted(glob.glob(root_folder+'*grid_V.nc'))


fig = plt.figure(figsize=(12,8))
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
ax = plt.axes(projection=ccrs.Mercator(central_longitude=180))                        
axpos = ax.get_position()                                                             
pos_x = axpos.x0+axpos.width + 0.01# + 0.25*axpos.width                               
pos_y = axpos.y0                                                                      
cax_width = 0.04                                                                      
cax_height = axpos.height                                                             
# Draw coastlines so we know where we are:                                            
ax.coastlines()                                                                       
# Set the map extent, making sure to specify the correct coordinate system            
# for geographical coordinates:                                                       
ax.set_extent([27.4, 42, 40.9, 47.5], crs=ccrs.PlateCarree())                          
ax.add_feature(cartopy.feature.LAND,color='gray')                                     

cbar_ax.set_position([axpos.x0 + axpos.width + 0.14, axpos.y0+0.04,
                          0.04, axpos.height*0.9])



# fig, ax1 = plt.subplots(figsize=(12,8)) 
ind = 0
for ind3, fnamet in enumerate(ls1t[60:80]):
    dft = xr.open_dataset(fnamet)
    fnameu = fnamet[:-4] + 'U.nc'
    fnamev = fnamet[:-4] + 'V.nc'
    dfu = xr.open_dataset(fnameu)
    dfv = xr.open_dataset(fnamev)
    zeta = rel_vort(dfu, dfv, gr)
    for ind2 in range(0,6):
        im = ax.pcolormesh(gr.nav_lon,gr.nav_lat,zeta[:,:,ind2].where(dft.tos[-1,:,:]!=0)/fcor,
                       cmap='BrBG_r',vmin=-0.6,vmax=0.6,transform=ccrs.PlateCarree());
        plt.colorbar(im,cax=cbar_ax)
        # plt.title(str(ind))
        plt.title(str(dfv.time_counter[ind2].coords)[-19:-9])
        print(ind)
        ind += 1
        fout = 'gifs/BS_rel_vort_' + str(ind).zfill(3) + '.png'          
        plt.savefig(fout, bbox_inches='tight',format='png',dpi=300) 
        # plt.clf()
    
