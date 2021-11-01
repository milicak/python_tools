import numpy as np
import subprocess
import matplotlib.pyplot as plt
from fem_tools import *
import sys

if len(sys.argv) < 3:
	print ('Error: Missing argument')
	print ('Please provide path of shyfem netcdf unstructured output ')
	print ('and path of partitioning file')
	sys.exit(1)		

shyfem_coord_path 	= sys.argv[1]
part_path 		= sys.argv[2]

lons,lats,elems,depth,time = load_shyfem_coords(shyfem_coord_path)


def read_partition(part_path):
	# read all lines apart from last
	data = np.genfromtxt(part_path,skip_header=1,skip_footer=1)
	# read last line in file
	line = subprocess.check_output(['tail', '-1', part_path])
	line = np.asarray([int(i) for i in line.split()])
	nr,nc = data.shape 
	data = data.reshape(nr*nc)
	# concatenate last line
	data = np.concatenate((data,line))
	return data

## read partition file
mpi_proc = read_partition(part_path)

mpi_proc        -= 1

nkn	= len(lons)
nel	= elems.shape[0]

nel_part	= len(mpi_proc)


if (nel != nel_part):
	print ('ERROR')	
	print ('Number of elements in shyfem grid does not match the')
	print ('number elements in partitioning file')
	sys.exit(1)


## this is to view process number at mouse 
x_centroid = np.mean(lons[elems],axis=1)
y_centroid = np.mean(lats[elems],axis=1)

glob_node	= np.arange(nkn)+1
glob_ele	= np.arange(nel)+1

def fmt(x, y):
        dist = np.sqrt((x-x_centroid)**2+(y-y_centroid)**2)
        dist2 = np.sqrt((x-lons)**2+(y-lats)**2)	
        idx = np.argmin(dist)
        z = mpi_proc[idx]
        e = glob_ele[idx]	
        idx2	= np.argmin(dist2)		
        n	= glob_node[idx2]	
        return 'x={x:.7f}  y={y:.7f} PROC={z:.1f} ele={e:6d} node={n:6d}'.format(x=x, y=y, z=z, e=e, n=n)

## make plot
plt.figure(figsize=(10,10))
cf = plt.tripcolor(lons,lats,elems,mpi_proc,cmap='hsv')
plt.triplot(lons,lats,elems,color='k',linewidth=0.2)
cb = plt.colorbar(cf)
cb.set_label('MPI Process',fontsize=15)
plt.gca().format_coord = fmt

plt.show()




