from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt

f = Dataset('ecmwf.nc','r')
precip 	= f.variables['precip'][:]
time 	= f.variables['time'][:]
lons 	= f.variables['lon'][:]
lats 	= f.variables['lat'][:]
f.close()

print 'precipitation dataset shape: ', precip.shape

subset = True

# box for subset
lon1, lon2 = 20, 32
lat1, lat2 = 36, 44

if subset == True:
	#latitude lower and upper index
	latli = np.argmin( np.abs( lats - lat1 ) )
	latui = np.argmin( np.abs( lats - lat2 ) ) 

	# longitude lower and upper index
	lonli = np.argmin( np.abs( lons - lon1 ) )
	lonui = np.argmin( np.abs( lons - lon2 ) )  

	# subset lons, lats, precip
	lons 	= lons[lonli:lonui]
	lats 	= lats[latui:latli]  # subsetting for lats is reversed
	precip 	= precip[:,latui:latli,lonli:lonui]

# invert along lat axis
precip = precip[:,::-1,:]

# set to 0 negative values
neg = np.where(precip < 0.0)
precip[neg] = 0.0

# convert rain from cumulated [m/6h]
# in mm/day 

# to convert to mm/day
conv = 1000.0 * 4
# not to convert and original
#conv = 1
precip = conv * precip

# convert time in seconds and create datetime list
time = [3600*i for i in time]
date0 = datetime(2017,1,1,0,0,0)
dates = [date0 + timedelta(seconds=i) for i in time]

nt,ny,nx = precip.shape
dx = abs(lons[1] - lons[0])
dy = abs(lats[1] - lats[0])
minlon = np.min(lons)
minlat = np.min(lats)
maxlon = np.max(lons)
maxlat = np.max(lats)

if subset == True:
	print 'subsetting to %s <= lon <= %s ' % (minlon,maxlon)
	print 'subsetting to %s <= lat <= %s ' % (minlat,maxlat)

date_format = '%Y%m%d %H%M%S\n'
header ='0 2 957839 %s 1 1 11\n' % (ny*nx)

## write to file
if conv==1:
    g = open('rain_utss_noconv.dat','w')
else:
    g = open('rain_utss_conv.dat','w')


for i in range(len(dates)):
	g.write(header)
	g.write(dates[i].strftime(date_format))
	g.write('%s %s %s %s %s %s 1e+20\n'%(nx,ny,minlon,minlat,dx,dy))
	g.write('rain [mm/day]\n')
	for y in range(ny):
		for x in range(nx):
			g.write('%f ' % precip[i,y,x])
			if x == nx-1:
				g.write('\n')

g.close()

#plt.imshow(precip[0,:])
#plt.show()

