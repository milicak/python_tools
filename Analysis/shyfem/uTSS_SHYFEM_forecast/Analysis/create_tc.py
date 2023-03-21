from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np


# root_folder = '/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_010/analysis/6h/netcdf/'
root_folder = '/data/inputs/metocean/rolling/atmos/ECMWF/IFS_010/1.0forecast/1h_3h_6h/netcdf/'
freq_ecmwf = 1
outfile = 'tc.dat'

date0 = datetime(2022,9,10,1,0,0)
subset = True
# box for subset
lon1, lon2 = 20, 32
lat1, lat2 = 36, 44

first_call = True

for itr in range(0,15): # 10 days
    date = date0 + timedelta(days=itr)
    ls1 = sorted(glob.glob(root_folder + date.strftime("%Y%m%d") + '/*MEDATL*fc00*'))
    df = xr.open_mfdataset(ls1)   
    lons = np.copy(df.lon)  
    lats = np.copy(df.lat) 

    if subset == True & first_call == True:
        #latitude lower and upper index
        latli = np.argmin( np.abs( lats - lat1 ) )
        latui = np.argmin( np.abs( lats - lat2 ) ) 
    
        # longitude lower and upper index
        lonli = np.argmin( np.abs( lons - lon1 ) )
        lonui = np.argmin( np.abs( lons - lon2 ) )  
    
    if subset == True:
        # subset lons, lats, precip
        lons    = lons[lonli:lonui]
        lats    = lats[latui:latli]  # subsetting for lats is reversed
        # only get the first day LATER we will change this
        airt = df.T2M[0:24,latui:latli,lonli:lonui]
        d2t     = df.D2M[0:24,latui:latli,lonli:lonui]
        tcc     = df.TCC[0:24,latui:latli,lonli:lonui]

    # invert along lat axis
    airt        = airt[:,::-1,:]
    d2t         = d2t[:,::-1,:]
    tcc         = tcc[:,::-1,:]
    airt = np.copy(airt)
    d2t  = np.copy(d2t)
    tcc  = np.copy(tcc)
    airt += -273.15
    d2t  += -273.15
    tcc  /= 100.0
    zerosolar = airt*0.0
    first_call = False
    
    # convert time in seconds and create datetime list
    if itr==0:
        time = [3600*freq_ecmwf*i for i in np.linspace(0,23,24)] 
        # time = [3600*freq_ecmwf*i for i in np.linspace(0,df.time.shape[0]-1,df.time.shape[0])] 
        dates = [date0 + timedelta(seconds=i) for i in time] 
    else:
        time = [(3600*freq_ecmwf*i+3600*freq_ecmwf) for i in np.linspace(0,23,24)]
        # time = [(3600*freq_ecmwf*i+3600*freq_ecmwf) for i in np.linspace(0,df.time.shape[0]-1,df.time.shape[0])]
        dates = [dates[-1] + timedelta(seconds=i) for i in time] 
    
    print(dates)
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
    header ='0 2 957839 %s 1 1 11\n' % (ny*nx)
    
    #     character*80, save :: srad = 'solar radiation [W/m**2]'
    #     character*80, save :: tair = 'air temperature [C]'
    #     character*80, save :: rhum = 'humidity [%]'
    #     character*80, save :: ccov = 'cloud cover [0-1]'
    #     character*80, save :: wbtm = 'wet bulb temperature [C]'
    #     character*80, save :: dewp = 'dew point temperature [C]'    
    # ## write to file
    g = open(outfile,'a')
    for i in range(len(dates)):
        g.write(header)
        g.write(dates[i].strftime(date_format))
        g.write('%s %s %s %s %s %s 1e+20\n'%(nx,ny,minlon,minlat,dx,dy))
        g.write('solar radiation [W/m**2]\n')
        for y in range(ny):
                for x in range(nx):
                        g.write('%f ' % zerosolar[i,y,x])
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
    
    
