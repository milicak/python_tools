import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,Point
from netCDF4 import Dataset

def distance_point_line(p1,p2,p):
	### returns the distance of point p from the line
	### represented by the 2 points p1,p2

	# p1 = (x1,y1) tuple with coordinates of 1st point on the line
	# p2 = (x2,y2) tuple with coordinates of 2st point on the line
	# p  = (xp,yp) tuple with coordinates of point to calculate distance

	x1 = float(p1[0])
	y1 = float(p1[1])
	x2 = float(p2[0])
	y2 = float(p2[1])
	xp = float(p[0])
	yp = float(p[1])

	if x1 == x2:
		# case of vertical line
		distance = np.abs(xp - x1)

	elif y1 == y2:
		# case of horizontal line
		distance = np.abs(yp - y1) 

	else:
		# normal case
		# determine line in the form a*x + b*y + c = 0

		a = 1. / (x2-x1)
		b = -1. / (y2-y1)
		c = -x1 / (x2-x1) + y1 / (y2-y1)

		distance = np.abs((a * xp + b * yp + c)) / np.sqrt(a**2 + b**2)

	#debug plot
	#fig,ax = plt.subplots()
	#ax.plot([x1,x2],[y1,y2])
	#ax.plot(xp,yp,'p')
	#ax.set_xlim([-5,5])
	#ax.set_ylim([-5,5])
	#ax.grid()
	#ax.set_aspect('equal')
	#plt.show()

	return distance	 	

def read_polygon_from_file(path):
        ## returns shapely.Polygon object from coordinates in file
        data = np.genfromtxt(path,usecols=(0,1),delimiter=',')
        x,y = data[:,0],data[:,1]
        poly = Polygon([xx,yy] for xx,yy in zip(x,y))
        return poly

def in_polygon(lons,lats,polygon):
        ### returns boolean array of (lons,lats) points inside/outside given shapely Polygon
        in_poly = []
        for x,y in zip(lons,lats):
                in_poly.append(polygon.contains(Point(x,y)))
        return np.asarray(in_poly)

def read_unstructured_dataset():
	## read SHYFEM
	return	

def interpolate_profile():
	## interpolates vertical profile to target grid 
	return

def read_regular_dataset(path,varname,target_depth):
	## read MFS or BSFS dataset 
	## and convert to scatter vertical profiles

	##  --------------------------------
	##  loading dataset
	##  ---------------------------------

	nc = Dataset(path)
	lon = nc.variables['lon'][:]
	lat = nc.variables['lat'][:]
	depth = nc.variables['depth'][:]
	variable = nc.variables[varname][:] 
	nc.close()

	##  --------------------------------
	##  converting to scatter of profiles
	##  ---------------------------------

	lon2,lat2 	= np.meshgrid(lon,lat)

	variable 	= variable.squeeze()

	nz,ny,nx 	= variable.shape
	
	## (nz,npoints2d)
	variable_profiles = variable.reshape(nz,nx*ny)

	# coordinares to 1D
	lat_1d 		= lat2.reshape(nx*ny)
	lon_1d 		= lon2.reshape(nx*ny)

	# 2D mask 
	mask2d 		= variable[0,:].mask
	mask1d 		= mask2d.reshape(nx*ny)
	inv_mask1d 	= mask1d == False

	# (nz,npoints2d) -> (npoints2d,nz)
	variable_profiles2 	= variable_profiles.T

	# select only sea profiles
	variable_profiles3 	= variable_profiles2[inv_mask1d]	
	lon_1d, lat_1d 		= lon_1d[inv_mask1d],lat_1d[inv_mask1d]

	n_profiles = len(lon_1d)

	##  --------------------------------
	##  Interpolating profiles onto target grid
	##  ---------------------------------
	
	variable_profiles4 	= np.zeros((n_profiles,len(target_depth)))

	for pp in range(n_profiles):
		good_profile 	= variable_profiles3[pp,variable_profiles3[pp,:].mask == False]
		good_depth 	= depth[variable_profiles3[pp,:].mask == False]
		variable_profiles4[pp,:] = np.interp(target_depth,good_depth,good_profile,right=good_profile[-1])

	return (lon_1d,lat_1d,variable_profiles4)

def load_profiles_dataset(path):
	nc = Dataset(path)
	lons = nc.variables['longitude'][:]
	lats = nc.variables['latitude'][:]
	temp = nc.variables['temperature'][:].squeeze()
	salt = nc.variables['salinity'][:].squeeze()
	depth = nc.variables['depth'][:]
	nc.close()
	return (lons,lats,depth,temp,salt)

def load_shyfem_dataset(path):
	nc = Dataset(path)
	lons = nc.variables['longitude'][:]
	lats = nc.variables['latitude'][:]
	temp = nc.variables['temperature'][:].squeeze()
	salt = nc.variables['salinity'][:].squeeze()
	depth = nc.variables['level'][:]
	nc.close()
	return (lons,lats,depth,temp,salt)

def interpolate_profiles_to_target_vertical_grid(profiles,original_depth,target_depth,degradation=1):
	profiles = profiles[::degradation,:]
	n_profiles = profiles.shape[0]
	# mask where profiles == 0
	profiles = np.ma.masked_where(profiles == 0 , profiles)
	profiles2 	= np.zeros((n_profiles,len(target_depth)))
	for pp in range(n_profiles):
		good_profile	= profiles[pp,profiles[pp,:].mask == False]
		good_depth 	= original_depth[profiles[pp,:].mask == False]
		#profiles2[pp,:] = np.interp(target_depth,good_depth,good_profile,right=good_profile[-1])	
		profiles2[pp,:] = np.interp(target_depth,good_depth,good_profile,right=0)

	profiles2 = np.ma.masked_where(profiles2 == 0, profiles2) 	
		
	return profiles2		
	

def extend_last_value(array):
	## it works only if you have this kind of situation
	## val val val val  0  0  0  0  0  0  0 
	## to get
	## val val val val val val val val val
	## Not for cases like
	## val  0 val val  0 val  0  0 

	#array2 = np.copy(array)
	n_profiles, n_levs = array.shape
	for pp in range(n_profiles):
		count = np.sum(array[pp] != 0)
		#count = np.ma.count(array[pp])
		array[pp,count:] = array[pp,count-1]

	return array
