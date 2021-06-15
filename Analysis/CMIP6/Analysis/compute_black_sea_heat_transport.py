import gcsfs
import geopy.distance

def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds

##################################################################################
###############

df = pd.read_csv('https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv')
dd = df[ (df.variable_id == 'thetao') & (df.experiment_id == 'ssp585') & (df.source_id == 'CNRM-ESM2-1')]
dt = read_data(dd.zstore.iloc[0])
dd = df[ (df.variable_id == 'vo') & (df.experiment_id == 'ssp585') & (df.source_id == 'CNRM-ESM2-1')]
dv = read_data(dd.zstore.iloc[0])
dz = np.copy(dt.lev_bounds[:,1])-np.copy(dt.lev_bounds[:,0])
Tsect = dt.thetao[:,:,210,315]
Vsect = dv.vo[:,:,210,316]

rho0=1025
Cp=3996

coords_1=(np.copy(dt.lat[210,315]),np.copy(dt.lon[210,315]))
coords_2=(np.copy(dt.lat[210,316]),np.copy(dt.lon[210,316]))
dy = geopy.distance.geodesic(coords_1, coords_2).m
# HT = int(rho_0*Cp * Tsect*Vsect*dz*dy)
HT  = rho0*Cp*Tsect*Vsect*dz*dy
HTsum = HT.sum('lev')
