import numpy as np
from netCDF4 import Dataset

def load_shyfem_coords(path):
	nc 	= Dataset(path)
	lons  	= nc.variables['longitude'][:]
	lats  	= nc.variables['latitude'][:]
	elems 	= nc.variables['element_index'][:]-1
	levels	= nc.variables['level'][:]
	time	= nc.variables['time'][:]
	nc.close()

	return (lons,lats,elems,levels,time)	 
	
def load_shyfem_variable(path,varname):
	nc 	= Dataset(path)
	var  	= nc.variables[varname][:]
	nc.close()

	return var	 

