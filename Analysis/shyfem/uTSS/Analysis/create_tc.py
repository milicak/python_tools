from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np

dataset = '2017_ecmwf_reanalysis_025.nc'
dataset = '2017_fv10_medatl.nc'
dataset = '2017_ifs_0125_analysis_medatl.nc'

nvars	= 4

f = Dataset(dataset,'r')
airt 	= f.variables['T2M'][:]
d2t 	= f.variables['D2M'][:]
tcc	= f.variables['TCC'][:]
time 	= f.variables['time'][:]
lons 	= f.variables['lon'][:]
lats 	= f.variables['lat'][:]
f.close()

airt += -273.15

print ('ecmwf dataset shape: ', airt.shape)

subset = True
subset = False

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

	# subset lons, lats, variables
	lons 	= lons[lonli:lonui]
	lats 	= lats[latui:latli]  # subsetting for lats is reversed
	windx 	= windx[:,latui:latli,lonli:lonui]
	windy 	= windy[:,latui:latli,lonli:lonui]
	patm 	= patm[:,latui:latli,lonli:lonui]

# invert along lat axis
airt 	= airt[:,::-1,:]
d2t 	= d2t[:,::-1,:]
tcc 	= tcc[:,::-1,:]


# convert time in seconds and create datetime list
time = [3600*i for i in time]
date0 = datetime(2017,1,1,0,0,0)
dates = [date0 + timedelta(seconds=i) for i in time]

nt,ny,nx = airt.shape
dx = abs(lons[1] - lons[0])
dy = abs(lats[1] - lats[0])
minlon = np.min(lons)
minlat = np.min(lats)
maxlon = np.max(lons)
maxlat = np.max(lats)

if subset == True:
	print ('subsetting to %s <= lon <= %s ' % (minlon,maxlon))
	print ('subsetting to %s <= lat <= %s ' % (minlat,maxlat))

date_format = '%Y%m%d %H%M%S\n'
header ='0 2 957839 %s 1 %d 11\n' % ((ny*nx),nvars)

## write to file
g = open('tc_SES_2017.dat','w')
for i in range(len(dates)):
	g.write(header)
	g.write(dates[i].strftime(date_format))
	g.write('%s %s %s %s %s %s 1e+20\n'%(nx,ny,minlon,minlat,dx,dy))
	g.write('solar radiation [W/m**2]\n')
	for y in range(ny):
		for x in range(nx):
			g.write('%f ' % airt[i,y,x])
			if x == nx-1:
				g.write('\n')
	g.write('air temperature [C]\n')
	for y in range(ny):
		for x in range(nx):
			g.write('%f ' % airt[i,y,x])
			if x == nx-1:
				g.write('\n')
	g.write('dew point temperature [C]\n')
	for y in range(ny):
		for x in range(nx):
			g.write('%f ' % d2t[i,y,x])
			if x == nx-1:
				g.write('\n')
	g.write('cloud cover [0-1]\n')
	for y in range(ny):
		for x in range(nx):
			g.write('%f ' % tcc[i,y,x])
			if x == nx-1:
				g.write('\n')

g.close()


