import utm
import scipy.io
import geopy.distance
import xmitgcm

time = 100

ds = xmitgcm.open_mdsdataset('/okyanus/users/milicak/models/MITgcm/Projects/urmia/work/',
                             prefix='SALT',geometry='cartesian',read_grid='True',
                             ref_date='2018-12-31 12:0:0',delta_t=2)

df = pd.read_excel('/okyanus/users/milicak/dataset/urmia/Water_density_section_April.xls')

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
    print(coords_2)
    dst[ind+1] = geopy.distance.geodesic(coords_1, coords_2).km



Ssurf = np.zeros(lon.shape)
for ind in range(0,lon.shape[0]):
    Ssurf[ind] = ds.SALT[time,0,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')

Ssect = np.zeros((ds.SALT.shape[1],lon.shape[0]))
for ind in range(0,lon.shape[0]):
    Ssect[:,ind] = ds.SALT[time,:,:,:].interp(XC=lon[ind],YC=lat[ind], method='linear')

Sbottom = np.array([Ssect[12,0],Ssect[16,1],
                    Ssect[17,2],Ssect[18,3],Ssect[15,4],Ssect[13,5],Ssect[15,6],
                    Ssect[14,7],Ssect[15,8],Ssect[14,9],Ssect[14,10],Ssect[13,11],
                    Ssect[14,12],Ssect[14,13],Ssect[14,14],Ssect[15,15],
                    Ssect[16,16],Ssect[16,17],Ssect[15,18],Ssect[15,19] ])


# compute mitgcm EOS using salinity and off set of 20 psu
rho_surf_mitgcm = 1237*(9.5e-4*(Ssurf+20-278))+(1237-1000)
rho_bottom_mitgcm = 1237*(9.5e-4*(Sbottom+20-278))+(1237-1000)

