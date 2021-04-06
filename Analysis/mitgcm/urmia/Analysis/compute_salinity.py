import utm
import scipy.io

df = pd.read_excel('/okyanus/users/milicak/dataset/urmia/Shallow_water_density_April_2019.xlsx')

lat, lon = utm.to_latlon(df.POINT_X, df.POINT_Y, 38, 'S')

lat = np.copy(lat)
lon = np.copy(lon)

# aa = df.columns.tolist()
# dens = np.copy(df[aa[-1]])
dens = np.copy(df.Value)

tt = np.concatenate((lon,lat,dens))
tt = np.reshape(tt,(3,1048575))
dens_april_shallow = tt
scipy.io.savemat('dens_april_shallow.mat', {'dens_april_shallow': dens_april_shallow})

df = pd.read_excel('/okyanus/users/milicak/dataset/urmia/Deep_water_density_April_2019.xlsx')

lat, lon = utm.to_latlon(df.POINT_X, df.POINT_Y, 38, 'S')

lat = np.copy(lat)
lon = np.copy(lon)

# aa = df.columns.tolist()
# dens = np.copy(df[aa[-1]])
dens = np.copy(df.value)

tt = np.concatenate((lon,lat,dens))
tt = np.reshape(tt,(3,1048575))
dens_april_deep = tt
scipy.io.savemat('dens_april_deep.mat', {'dens_april_deep': dens_april_deep})
