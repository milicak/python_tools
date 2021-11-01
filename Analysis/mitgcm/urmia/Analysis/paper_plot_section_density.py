import utm
import scipy.io
import geopy.distance
import xmitgcm


def find_rho(ds,dt):
    rhoNil = 1237
    sBeta = 9.5e-4
    tAlpha = 3.5e-4
    dRho = 1237 - 1000
    refSalt = 278
    refTemp = 23
    rho = rhoNil*(sBeta*(ds.SALT-refSalt)-tAlpha*(dt.THETA-refTemp) ) + dRho
    return rho

def find_rho2(ds):
    rhoNil = 1237
    sBeta = 9.5e-4
    tAlpha = 3.5e-4
    dRho = 1237 - 1000
    refSalt = 278
    refTemp = 23
    rho = rhoNil*(sBeta*(ds.SALT-refSalt)-tAlpha*(18-refTemp) ) + dRho
    return rho

ds = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/', prefix='SALT',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work/', prefix='THETA',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
ds2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/', prefix='SALT',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt2 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work2/', prefix='THETA',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
ds3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/', prefix='SALT',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')
dt3 = xmitgcm.open_mdsdataset('/work/opa/mi19918/Projects/mitgcm/urmia/work3/', prefix='THETA',
                        ref_date = "2018-12-31 0:0:0", delta_t=2, geometry='cartesian', read_grid='True')

rho = find_rho(ds,dt)
rho2 = find_rho(ds2,dt2)
rho3 = find_rho(ds3,dt3)

rho = find_rho2(ds)
rho2 = find_rho2(ds2)
rho3 = find_rho2(ds3)

rho = rho.where(rho>0,0)
rho2 = rho2.where(rho2>0,0) 
rho3 = rho3.where(rho3>0,0) 

df = pd.read_excel('Water_density_section_April.xls')

lat, lon = utm.to_latlon(df.Lat, df.Lon, 38, 'S')

rho_deep = df['Density of \ndeep water']
rho_surface = df['Density of \nShallow water']

lat = np.copy(lat)
lon = np.copy(lon)

dst = np.zeros(lon.shape)                                       
coords_1 = (lat[0], lon[0])                                     
for ind in range(0,lon.shape[0]-1):                             
    # coordinates of M23 (lat, lon)                             
    coords_2 = (lat[ind+1], lon[ind+1])                         
    dst[ind+1] = geopy.distance.geodesic(coords_1, coords_2).km 
                                                                

time = 110 #99 
rhosurf = np.zeros(lon.shape)                                                     
rhosurf2 = np.zeros(lon.shape)                                                     
rhosurf3 = np.zeros(lon.shape)                                                     
for ind in range(0,lon.shape[0]):                                               
    rhosurf[ind] = rho[time,0,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')
    rhosurf2[ind] = rho2[time,0,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')
    rhosurf3[ind] = rho3[time,0,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')


Ssect = np.zeros((ds.SALT.shape[1],lon.shape[0]))                               
Ssect2 = np.zeros((ds.SALT.shape[1],lon.shape[0]))                               
Ssect3 = np.zeros((ds.SALT.shape[1],lon.shape[0]))                               
for ind in range(0,lon.shape[0]):                                               
    Ssect[:,ind] = rho[time,:,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')
    Ssect2[:,ind] = rho2[time,:,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')
    Ssect3[:,ind] = rho3[time,:,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')
                                                                                
rhobottom = np.array([Ssect[12,0],Ssect[16,1],                                    
                    Ssect[17,2],Ssect[18,3],Ssect[15,4],Ssect[13,5],Ssect[15,6],
                    Ssect[14,7],Ssect[15,8],Ssect[14,9],Ssect[14,10],Ssect[13,11],
                    Ssect[14,12],Ssect[14,13],Ssect[14,14],Ssect[15,15],        
                    Ssect[16,16],Ssect[16,17],Ssect[15,18],Ssect[15,19] ])      
                                                                                
rhobottom2 = np.array([Ssect2[12,0],Ssect2[16,1],                                    
                    Ssect2[17,2],Ssect2[18,3],Ssect2[15,4],Ssect2[13,5],Ssect2[15,6],
                    Ssect2[14,7],Ssect2[15,8],Ssect2[14,9],Ssect2[14,10],Ssect2[13,11],
                    Ssect2[14,12],Ssect2[14,13],Ssect2[14,14],Ssect2[15,15],        
                    Ssect2[16,16],Ssect2[16,17],Ssect2[15,18],Ssect2[15,19] ])      

fig, axes = plt.subplots(figsize=(10,6))
ax1 = plt.subplot2grid(shape=(1,1),loc=(0,0), colspan=1)
ax1.plot(dst,rho_surface-1e3,'-*',label='obs-surface-density')  
ax1.plot(dst,rho_deep-1e3,'-*',label='obs-bottom-density')  
ax1.plot(dst,rhosurf+20,'-*',label='run1-surface-density')  
ax1.plot(dst,rhosurf+39,'-*',label='run1-bottom-density')  
ax1.plot(dst,rhosurf2+20,'-*',label='run3-surface-density')  
ax1.plot(dst,rhosurf2+22,'-*',label='run3-bottom-density')  
ax1.legend()
ax1.set_xlabel('Distance [km]', fontsize=14);
ax1.set_ylabel('Density anomaly [$kgm^{-3}$]', fontsize=14);
x = list(range(0, 100, 20))
ax1.set_xticklabels(x, fontsize=14);
y = list(range(60, 230, 20))
ax1.set_yticklabels(y, fontsize=14);

plt.savefig('paperfigs/density_vertical_section.png', bbox_inches='tight',format='png',dpi=300)


