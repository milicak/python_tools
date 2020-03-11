import numpy as np
from init_utils import distance_point_line, read_regular_dataset
from init_utils import load_profiles_dataset, load_shyfem_dataset
from init_utils import interpolate_profiles_to_target_vertical_grid 
from init_utils import extend_last_value
from shapely.geometry import Point, Polygon, LineString
from grd_tools import read_grd_nodes
from scipy.spatial import cKDTree
from datetime import datetime
#import matplotlib.pyplot as plt
import sys

startTime = datetime.now()

print 'reading target horizontal grid '
### read target grid 
grd_path = 'ses_v6c_hmin3_hmax4500.grd'
node_idx,xgrid,ygrid = read_grd_nodes(grd_path,mode=0)

nkn = len(xgrid)

### read target vertical grid
shyfem_bottom_levels = [2, 4, 6, 8, 10, 12,  14,  16,  18,  20,  22,  24,  26,  28,  30,  32,  34,  36,  38,  40
         , 42,  44,  46,  48,  50,  52,  54,  56,  58,  60 , 62,  64, 66, 68,  70,  75, 80, 85
         , 90, 95, 100, 105, 110,  120,  130,  140,  150,  160,  170,  180,  190,  200
         , 220,  240,  260, 280,  300,  320,  340,  360,  380,  400,  420,  440,  460,  480,  500,  550
         ,  600,  650,  700,  750,  800,  850,  900,  950,  1000,  1100,  1200,  1300,  1400 ,  1500
         , 1600,  1700,  1800, 1900,  2000,  2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000
         , 3150,  3300, 3500, 3700, 3900, 4200, 4500]
edges = np.asarray([0]+shyfem_bottom_levels)
center_levels = 0.5*(edges[1:]+edges[:-1])

n_target_levs = len(center_levels)

###################################
### LOADING DATASETS TO INTERPOLATE
###################################	

print 'loading dataset / converting to profiles / interpolate onto target vertical grid ...'

### read MFS dataset
#mfs_path = '201601_mfs_profiles_interpolated.nc'
mfs_path = '2016_2017_mfs_profiles_interpolated.nc'
mfs_lons,mfs_lats,mfs_depth,mfs_temp,mfs_salt = load_profiles_dataset(mfs_path)
#bsfs_path = '201601_bsfs_profiles_interpolated.nc'
bsfs_path = '2016_2017_bsfs_profiles_interpolated.nc'
bsfs_lons,bsfs_lats,bsfs_depth,bsfs_temp,bsfs_salt = load_profiles_dataset(bsfs_path)

print mfs_temp.shape
print bsfs_temp.shape

### load uTSS dataset (yet to interpolate onto target vertical grid)
#utss_path = '201601_utss.nc'
utss_path = '2016_2017_utss.nc'
utss_lons, utss_lats, utss_depth, utss_temp,utss_salt = load_shyfem_dataset(utss_path)
degradation = 2
utss_temp = interpolate_profiles_to_target_vertical_grid(utss_temp,utss_depth,center_levels,degradation=degradation)
utss_salt = interpolate_profiles_to_target_vertical_grid(utss_salt,utss_depth,center_levels,degradation=degradation)
utss_lons = utss_lons[::degradation]
utss_lats = utss_lats[::degradation]


print utss_temp.shape

print 'creating KDTrees from datasets ..'
### create KDTree for MFS
mfs_tree = cKDTree(zip(mfs_lons,mfs_lats))
### create KDTree for uTSS
utss_tree = cKDTree(zip(utss_lons,utss_lats))
### create KDTree for BSFS
bsfs_tree = cKDTree(zip(bsfs_lons,bsfs_lats))

grade = 4

print 'making query on KDTrees, grade = %s ....' % grade
### make query of closest MFS obs to target grid 	--> d_mfs, inds_mfs
d_mfs, inds_mfs = mfs_tree.query(zip(xgrid,ygrid),k = grade)
### make query of closest uTSS obs to target grid 	--> d_utss, inds_utss
d_utss, inds_utss = utss_tree.query(zip(xgrid,ygrid),k = grade)
### make query of closest BSFS obs to target grid	--> d_bsfs, inds_bsfs
d_bsfs, inds_bsfs = bsfs_tree.query(zip(xgrid,ygrid),k = grade)


#####################################
### BUFFER ZONES - alpha coefficients
#####################################


### Define points of buffer zone1 (MFS-uTSS)
A1 = (25,38.8)
A2 = (25,39.8)
A3 = (25,50)
A4 = (25,30)

p1 = Point(A1)
p2 = Point(A2)
p3 = Point(A3)
p4 = Point(A4)

### Define points of buffer zone1 (uTSS-BSFS)
A5 = (29.088,41.3025)
A6 = (29.337,41.916)
A7 = (32.0,48.46)
A8 = (23.68,28.01)


p5 = Point(A5)
p6 = Point(A6)
p7 = Point(A7)
p8 = Point(A8)

# distances between lines of 2 buffer areas
dist_p1_p2 = p1.distance(p2)
dist_p5_p6 = p5.distance(p6)

## bounding lines of buffer areas (2 for each buffer zone)
## first 2 point are on the 2 lines. The third is sufficiently "far".
## Note that the order of points defines the orientation of line

## from p1 to p2
line1 = LineString([p1,p2,p3])
## from p2 to p1
line2 = LineString([p2,p1,p4])
## from p5 to p6
line3 = LineString([p5,p6,p7])
## from p6 to p5
line4 = LineString([p6,p5,p8])

## initialize oriented distance of mesh points 
## projected on the lines

d1 = np.zeros(nkn)
d2 = np.zeros(nkn)
d3 = np.zeros(nkn)
d4 = np.zeros(nkn)

print 'Calculating distances along buffer zones..'
for nn in range(nkn):
	d1[nn] = line1.project(Point(xgrid[nn],ygrid[nn]))
	d2[nn] = line2.project(Point(xgrid[nn],ygrid[nn]))
	d3[nn] = line3.project(Point(xgrid[nn],ygrid[nn]))
	d4[nn] = line4.project(Point(xgrid[nn],ygrid[nn]))
	

alpha1 = -d1 / dist_p1_p2 + 1
alpha1[ alpha1 < 0] = 0
alpha2 = -d2 / dist_p1_p2 + 1
alpha2[ alpha2 < 0] = 0
alpha3 = -d3 / dist_p5_p6 + 1
alpha3[ alpha3 < 0] = 0
alpha4 = -d4 / dist_p5_p6 + 1
alpha4[ alpha4 < 0] = 0

alpha_mfs = alpha1
## correct alpha_mfs west of greece, it must be 1
alpha_mfs[xgrid < 22] = 1 
alpha_utss = alpha2 * alpha3 
## correct alpha_utss west of greece, it must be 0
alpha_utss[xgrid < 22] = 0
## correct alpha bsfs below turkey (there is a small area non zero)
alpha4[ygrid < 37] = 0
alpha_bsfs = alpha4

## check whether alphas are normalized
sum_alphas = np.sum( alpha_mfs + alpha_utss + alpha_bsfs)
sum_alphas = sum_alphas / nkn

print 'sum alphas ', sum_alphas
if sum_alphas != 1:
	print 'WARNING: the sum of alphas should be 1'

## expand alphas along neighbours
#alpha_mfs = np.expand_dims(alpha_mfs,axis=1)
#alpha_mfs = np.repeat(alpha_mfs,grade,axis=1)
#alpha_utss = np.expand_dims(alpha_utss,axis=1)
#alpha_utss = np.repeat(alpha_utss,grade,axis=1)
#alpha_bsfs = np.expand_dims(alpha_bsfs,axis=1)
#alpha_bsfs = np.repeat(alpha_bsfs,grade,axis=1)

## expand aplhas along depth
alpha_mfs = np.expand_dims(alpha_mfs,axis=1)
alpha_mfs = np.repeat(alpha_mfs,n_target_levs,axis=1)
alpha_utss = np.expand_dims(alpha_utss,axis=1)
alpha_utss = np.repeat(alpha_utss,n_target_levs,axis=1)
alpha_bsfs = np.expand_dims(alpha_bsfs,axis=1)
alpha_bsfs = np.repeat(alpha_bsfs,n_target_levs,axis=1)

print 'shape of alphas'
print alpha_mfs.shape
print alpha_utss.shape
print alpha_bsfs.shape

#np.savez('alphas.npz',alpha_mfs=alpha_mfs, alpha_utss=alpha_utss,alpha_bsfs=alpha_bsfs)


print 'Calculating weights based on distance..'

L = 0.25 # correlation distance

### calculate weights from d_mfs, d_utss, d_bsfs
w_mfs = np.exp(-(d_mfs/L)**2)
w_utss = np.exp(-(d_utss/L)**2)
w_bsfs = np.exp(-(d_bsfs/L)**2)

#sys.exit(0)


print 'first values of w_mfs'
print w_mfs[:10]


## I have to expand weights array to make
## multiplication. Basically I have to repeat the
## weights array along depth direction making this
## transformation of shape
## (nodes,neighbours) --> (nodes,neighbours,depth)

w_mfs = np.expand_dims(w_mfs,axis=2)
w_mfs = np.repeat(w_mfs,n_target_levs,axis=2)
w_utss = np.expand_dims(w_utss,axis=2)
w_utss = np.repeat(w_utss,n_target_levs,axis=2)
w_bsfs = np.expand_dims(w_bsfs,axis=2)
w_bsfs = np.repeat(w_bsfs,n_target_levs,axis=2)

# give a minimum value to weights
w_mfs = np.clip(w_mfs,a_min=0.01,a_max=1)
w_utss = np.clip(w_utss,a_min=0.01,a_max=1)
w_bsfs = np.clip(w_bsfs,a_min=0.01,a_max=1)

w_mfs  = np.ma.masked_where( mfs_temp[inds_mfs,:].mask == True  , w_mfs)
w_utss = np.ma.masked_where( utss_temp[inds_utss,:].mask == True, w_utss)
w_bsfs = np.ma.masked_where( bsfs_temp[inds_bsfs,:].mask == True, w_bsfs)

######################################################
### CALCULATE WEIGHTED AVERAGE for T,S
######################################################

#np.savez('temps.npz',mfs_temp=mfs_temp[inds_mfs].data,utss_temp=utss_temp[inds_utss].data,bsfs_temp=bsfs_temp[inds_bsfs].data)


print 'Calculating weighted average..'
mfs_part_temp  = np.sum( w_mfs  * mfs_temp[inds_mfs]  , axis = 1) / np.sum( w_mfs , axis = 1)
utss_part_temp = np.sum( w_utss * utss_temp[inds_utss], axis = 1) / np.sum( w_utss, axis = 1)
bsfs_part_temp = np.sum( w_bsfs * bsfs_temp[inds_bsfs], axis = 1) / np.sum( w_bsfs, axis = 1)
mfs_part_temp[mfs_part_temp.mask == True] = 0
utss_part_temp[utss_part_temp.mask == True] = 0
bsfs_part_temp[bsfs_part_temp.mask == True] = 0

mfs_part_salt  = np.sum( w_mfs  * mfs_salt[inds_mfs]  , axis = 1) / np.sum( w_mfs , axis = 1)
utss_part_salt = np.sum( w_utss * utss_salt[inds_utss], axis = 1) / np.sum( w_utss, axis = 1)
bsfs_part_salt = np.sum( w_bsfs * bsfs_salt[inds_bsfs], axis = 1) / np.sum( w_bsfs, axis = 1)
mfs_part_salt[mfs_part_salt.mask == True] = 0
utss_part_salt[utss_part_salt.mask == True] = 0
bsfs_part_salt[bsfs_part_salt.mask == True] = 0
# I think is time to extend the last good value to the whole water column


print 'extending last good value to whole column...'

mfs_part_temp2  = extend_last_value(mfs_part_temp) 
utss_part_temp2 = extend_last_value(utss_part_temp) 
bsfs_part_temp2 = extend_last_value(bsfs_part_temp) 

mfs_part_salt2  = extend_last_value(mfs_part_salt) 
utss_part_salt2 = extend_last_value(utss_part_salt) 
bsfs_part_salt2 = extend_last_value(bsfs_part_salt) 


temp_init = alpha_mfs * mfs_part_temp2 + alpha_utss * utss_part_temp2 + alpha_bsfs * bsfs_part_temp2 
salt_init = alpha_mfs * mfs_part_salt2 + alpha_utss * utss_part_salt2 + alpha_bsfs * bsfs_part_salt2 

#np.savez('temp_init_before_rearranging.npz',temp_init=temp_init,xgrid=xgrid,ygrid=ygrid)

#sys.exit(0)


######################################################
### REARRANGE FOR INTERNAL INTERNAL-EXTERNAL NUMBERING
#####################################################
print 're-arranging with internal ordering..'

in_ext_path = 'in_ext_nodes_ses_v6c.dat'
internal_idx = np.genfromtxt(in_ext_path,usecols=1,dtype=int)

rank = np.argsort(node_idx)
inds_reorder = rank[np.searchsorted(node_idx,internal_idx,sorter=rank)]
temp_init = temp_init[inds_reorder,:]
salt_init = salt_init[inds_reorder,:]
xgrid = xgrid[inds_reorder]
ygrid = ygrid[inds_reorder]

#print 'replacing Nones with last profile values'
#for nn in range(nkn):
#        ## get last non masked value
#        last_val = temp_init[nn,temp_init[nn,:].mask == False][-1]
#        #last_val_salt = salt_init[nn,salt_init[nn,:].mask == False][-1]
#        ## give this value to the masked values
#        temp_init[nn,temp_init[nn,:].mask == True] = last_val
#        #salt_init[nn,salt_init[nn,:].mask == True] = last_val_salt

###########################
### WRITE TO FILE
###########################


print 'writing files...'
### write to file
header = '0 2 957839 %d %d 1 1\n'
date_init = datetime(2016,1,1,0,0,0)
date_init_str = date_init.strftime('%Y%m%d %H%M%S\n')
bottom_levels_list = ' '.join(str(ii) for ii in shyfem_bottom_levels) + '\n'
varname = 'temperature [C]\n'
varname_salt = 'salinity [psu]\n'

#fout_temp = 'tempin.dat' 
#fout_salt = 'saltin.dat' 
fout_temp = 'tempin_2016_2017a.dat' 
fout_salt = 'saltin_2016_2017a.dat' 

ff = open(fout_temp,'w')
ff2 = open(fout_salt,'w')

ff.write(header % (nkn, n_target_levs))
ff.write(date_init_str)
ff.write(bottom_levels_list)
ff.write(varname)
ff2.write(header % (nkn, n_target_levs))
ff2.write(date_init_str)
ff2.write(bottom_levels_list)
ff2.write(varname_salt)
string0 = '%s -999.0 ' % n_target_levs
for nn in range(nkn):
        stringg = string0 + ' '.join(' %.5f ' % ii for ii in temp_init[nn,:].tolist()) + '\n'
        stringg_salt = string0 + ' '.join(' %.5f' % ii for ii in salt_init[nn,:].tolist()) + '\n'
        ff.write(stringg)
        ff2.write(stringg_salt)
ff.close()
ff2.close()


print 'total execution time: ', datetime.now() - startTime




