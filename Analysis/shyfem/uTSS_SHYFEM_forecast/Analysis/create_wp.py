from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np


# root_folder = '/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_010/analysis/6h/netcdf/'
root_folder = '/data/inputs/metocean/rolling/atmos/ECMWF/IFS_010/1.0forecast/1h_3h_6h/netcdf/'
freq_ecmwf = 1
outfile = 'wp.dat'

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
        windx   = df.U10M[0:24,latui:latli,lonli:lonui]
        windy   = df.V10M[0:24,latui:latli,lonli:lonui]
        patm    = df.MSL[0:24,latui:latli,lonli:lonui]

    # invert along lat axis
    windx = windx[:,::-1,:]
    windy = windy[:,::-1,:]
    patm  = patm[:,::-1,:]
    windx = np.copy(windx)
    windy = np.copy(windy)
    patm  = np.copy(patm)
    first_call = False
    
    # No need to convert to Pa because ECMWF has Pa in these years
    # convert patm to Pa
    # patm      = 100 * patm
    
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
    nt,ny,nx = windx.shape
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
    
    ## write to file
    g = open(outfile,'a')
    for i in range(len(dates)):
        g.write(header)
        g.write(dates[i].strftime(date_format))
        g.write('%s %s %s %s %s %s 1e+20\n'%(nx,ny,minlon,minlat,dx,dy))
        g.write('wind velocity in x [m/s]\n')
        for y in range(ny):
                for x in range(nx):
                        g.write('%f ' % windx[i,y,x])
                        if x == nx-1:
                                g.write('\n')
        g.write('wind velocity in y [m/s]\n')
        for y in range(ny):
                for x in range(nx):
                        g.write('%f ' % windy[i,y,x])
                        if x == nx-1:
                                g.write('\n')
        g.write('pressure (atmospheric) [Pa]\n')
        for y in range(ny):
                for x in range(nx):
                        g.write('%f ' % patm[i,y,x])
                        if x == nx-1:
                                g.write('\n')
    
    g.close()
    
