''' This subroutine computes the vertical section of a scalar (temperature,
salinity) along the pathway of Atlantic Water in the Nordic Seas. The section
location can be found in Ilicak et al. (2016) paper '''
import numpy as np
import sys
from netcdf_functions import nc_read

sys.path.insert(0,'/home/mil021/anaconda2/envs/esmpy/lib/python2.7/site-packages/')


# lon and lat values of the AW section
lon_AW = [17.6, 16.5, 16.05, 15.6, 15.1, 14.1, 13.0, 12.0, 10.0, 8.0, 4.0, 4.0,
        10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0,
        120.0, 130.0, 140.0]
lat_AW = [69.0, 70.6, 71.3, 72.02, 72.8, 73.8, 75.0, 76.0, 77.0, 78.0, 79.0,
          80.0, 81.0, 81.8, 81.8, 82.6, 83.0, 83.2, 83.1, 82.8, 82.5, 81.8,
          79.7, 78.2, 78.2, 79.7]
lon_AW = np.array(lon_AW)
lat_AW = np.array(lat_AW)


#section_name=[]
#section_name['Atlantic inflow']

gridfile='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
lonvar = 'plon'
latvar = 'plat'
loncvar = 'pclon'
latcvar = 'pclat'
areavar = 'parea'
maskvar = 'pmask'
sclfa = 1.0 # area scale factor 1 for m2 and 0.01 for cm2

# first you need to create grid file make_remap_grid
def make_remap_grid(gridfile, lonvar, latvar, loncvar, latcvar, areavar, maskvar, sclfa):
    ''' This function creates remapping grid files '''
    deg2rad = np.pi/180.0
    grid_center_lon = nc_read(gridfile, lonvar)
    grid_center_lat = nc_read(gridfile, latvar)
    grid_corner_lon = nc_read(gridfile, loncvar)
    grid_corner_lat = nc_read(gridfile, latcvar)
    grid_area = nc_read(gridfile, areavar)
    m2rad = 1.569612305760477e-07
    #convert m2 area to rad2
    grid_area = grid_area*m2rad*m2rad*sclfa
    ny = np.shape(grid_area)[0]
    nx = np.shape(grid_area)[1]
    grid_mask = nc_read(gridfile, maskvar)
    # map variables
    mapvars = {}
    mapvars['lon'] = grid_center_lon
    mapvars['lat'] = grid_center_lat
    mapvars['lonc'] = grid_corner_lon
    mapvars['latc'] = grid_corner_lat
    mapvars['area'] = grid_area
    mapvars['mask'] = grid_mask

    return mapvars



# call make_remap_grid
mapvars = make_remap_grid(gridfile, lonvar, latvar,
                                              loncvar, latcvar, areavar,
                                              maskvar, sclfa)


# second you need to create map_file using create_section_map
def create_section_map(grid_file, section_name, lonin, latin, mapvars):
    ''' This function creates section map '''




# third you need to create pcolor variables using map_scalar2section
