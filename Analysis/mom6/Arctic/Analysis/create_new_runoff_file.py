import netCDF4
from matplotlib import path

def crmask(lonlat_polygon, lonlat_model, lonmodel):
    p = path.Path(lonlat_polygon)
    mask = p.contains_points(lonlat_model)
    mask.shape = lonmodel.shape
    mask = np.multiply(mask, 1)
    return mask


lons_lats_poly= np.genfromtxt('Arctic_region_map.txt',usecols=(0,1))

df = xr.open_dataset('/okyanus/users/milicak/dataset/CORE2/NYF_v2.0/runoff.daitren.clim.v2011.02.10.nc')
ni = np.copy(df.ni)
nj = np.copy(df.nj)
lon, lat = np.meshgrid(ni, nj)

m = Basemap(projection='npstere',boundinglat=-89.5,lon_0=0,resolution='l')
x, y = m(lon,lat)
x1, y1 = m(lonarc,latarc)
lons_lats_poly_rot = np.column_stack((x1.flatten(),y1.flatten()))
lons_lats_model_rot = np.column_stack((x.flatten(),x.flatten()))
mask = crmask(lons_lats_poly_rot,lons_lats_model_rot,x)

lons_lats_model = np.column_stack((lon.flatten(),lat.flatten()))
mask = crmask(lons_lats_poly,lons_lats_model,lon)

dset = netCDF4.Dataset('runoff_daitren_clim_Arctic.nc','r+')

dset['area'][:]=area

dset.close()
