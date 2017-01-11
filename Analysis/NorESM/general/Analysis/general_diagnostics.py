#!/usr/bin/env python 
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

plt.ion()

datesep = '-'

# tripolar 1degree grid
grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
# tripolar 0.25degree grid
#grid_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/maps/grid_0_25degree.nc';
# bi-polar grid
#grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';

# global constants
aradius = 6.373E6      # Radius of Earth (m)
Lv = 2.5e6  # Joule/kg
rhow = 1000  # kg/m3
Lhvap = 2.5E6    # Latent heat of vaporization (J / kg)
Lhsub = 2.834E6   # Latent heat of sublimation (J / kg)
Lhfus = Lhsub - Lhvap  # Latent heat of fusion (J / kg)


def ncread_time_surface(fname, variable, timestr, timeend, x, y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp = np.zeros([y, x])
    for i in range(timestr, timeend):
        # print i
        tmp = tmp+ncfile.variables[variable][i, :, :].copy()

    tmp = tmp/(timeend-timestr)
    return tmp


def get_sdate_ini(root_folder,cmpnt,mdl,ext):
    # Get dimensions and time attributes for ocn or atm or others ...
    prefix=root_folder+expid+'/'+cmpnt+'/hist/'+expid+ '.'+ mdl + \
            '.'+ext+'.'
    if m2y==1:
        sdate="%4.4d%c%2.2d" % (fyear,datesep,1)
    else:
        sdate="%4.4d" % (fyear)


    return prefix, sdate


def inferred_heat_transport(energy_in, lat_deg):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate
    lat_rad = np.deg2rad(lat_deg)
    return (1E-15 * 2 * np.math.pi * aradius**2 *
            integrate.cumtrapz( np.cos(lat_rad)*energy_in, \
            x=lat_rad, initial=0.))


def timemean(prefix, var, varname, fyear, lyear):
    ''' timemean averages '''
    months2days=[31,  28,  31,  30,  31,   30,   31,  31,   30, 31,   30, 31];
    yeardays=sum(months2days);
    mw = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float)
    mw = mw/sum(mw)
    n=0.0;
    for year in xrange(fyear,lyear+1):
        if m2y==1:
            for month in xrange(0,12):
                n=n+mw[month]
                sdate="%4.4d%c%2.2d" % (year,datesep,month+1)
                dnm=nc_read(prefix+sdate+'.nc',varname);
                var=var+np.squeeze(dnm.data)*mw[month]
            print sdate,n
        
            
        else:
            sdate = "%4.4d" % (year)
            n += 1.0
            dnm=nc_read(prefix+sdate+'.nc',varname);
            var=var+np.squeeze(dnm.data)
            print sdate,n

    var = var/n
    return var



def timemean_old(args):
    ''' timemean averages '''
    mask=nc_read(grid_file,'pmask');

    prefix,sdate = get_sdate_ini(root_folder) 
    months2days=[31,  28,  31,  30,  31,   30,   31,  31,   30, 31,   30, 31];
    yeardays=sum(months2days);

    nx=ncgetdim(prefix+sdate+'.nc','x');
    ny=ncgetdim(prefix+sdate+'.nc','y');
    nz=ncgetdim(prefix+sdate+'.nc','depth');
    depth=nc_read(prefix+sdate+'.nc','depth');

    templvl=np.zeros((nz,ny,nx));
    salnlvl=np.zeros((nz,ny,nx));
    difisolvl=np.zeros((nz,ny,nx));
    difdialvl=np.zeros((nz,ny,nx));
    mask=np.tile(mask,(nz,1,1));
    #mask=np.transpose(mask, (2, 1, 0));

    n=0;
    if m2y==1:
         for year in xrange(fyear,lyear+1):
              for month in xrange(0,12):
                   n=n+months2days[month]
                   sdate="%4.4d%c%2.2d" % (year,datesep,month+1)
                   dnm=nc_read(prefix+sdate+'.nc','templvl');
                   templvl=templvl+np.squeeze(dnm.data)*mask*months2days[month]
                   dnm=nc_read(prefix+sdate+'.nc','salnlvl');
                   salnlvl=salnlvl+np.squeeze(dnm.data)*mask*months2days[month]
                   dnm=nc_read(prefix+sdate+'.nc','difisolvl');
                   difisolvl=difisolvl+10**(np.squeeze(dnm.data))*mask*months2days[month]
                   dnm=nc_read(prefix+sdate+'.nc','difdialvl');
                   difdialvl=difdialvl+10**(np.squeeze(dnm.data))*mask*months2days[month]
                   print sdate



def passagevolumetransporttime(root_folder,sectionno):
    ''' compute Drake Passage transport in time'''

    vol = {}
    sectiondr = 5 # Drake passage
    for year in xrange(fyear,lyear+1):
        if m2y==1:
            for month in xrange(0,12):
               sdate = "%4.4d%c%2.2d" % (year,datesep,month+1)
               dnm = nc_read(prefix+sdate+'.nc','voltr');
               secname = nc_read(prefix+sdate+'.nc','section');
               dnm = np.copy(dnm[:,sectionno-1])
               vol.setdefault('Volume',[]).append(dnm*1e-9)
               print sdate


        else:
            sdate = "%4.4d" % (year)
            dnm = nc_read(prefix+sdate+'.nc','voltr');
            dnm = np.copy(dnm[:,sectionno-1])
            secname = nc_read(prefix+sdate+'.nc','section');
            vol.setdefault('Volume',[]).append(dnm*1e-9)
            print sdate


    aa = ''.join(secname[sectionno-1])
    plt.plot(vol['Volume'],'b',label=aa)
    plt.legend(loc='lower right')
    return vol


def amoctime(root_folder,cmpnt,mdl,ext):
    ''' compute max amoc, amoc26 in time'''
    
    prefix,sdate = get_sdate_ini(root_folder,cmpnt,mdl,ext) 
    # Latitude interval for maximum AMOC search
    lat=nc_read(prefix+sdate+'.nc','lat');
    lat1=20;
    lat2=60;
    lat3=26.5;
    ind1=np.min(np.where(lat>=lat1));
    ind2=np.max(np.where(lat<=lat2));
    ind3=np.max(np.where(lat<=lat3));
    amoc = {}

    for year in xrange(fyear,lyear+1):
        if m2y==1:
            for month in xrange(0,12):
#              n=n+months2days[month]
               sdate = "%4.4d%c%2.2d" % (year,datesep,month+1)
               dnm = nc_read(prefix+sdate+'.nc','mmflxd');
               dnm = np.copy(dnm[0,region-1,:,:])
               dnm[dnm>1e25] = np.nan
               amoc.setdefault('max',[]).append(np.nanmax(dnm[:,ind1-1:ind2-1])*1e-9)
               amoc.setdefault('26N',[]).append(np.nanmax(dnm[:,ind3-1])*1e-9)
               print sdate


        else:
            sdate = "%4.4d" % (year)
            dnm = nc_read(prefix+sdate+'.nc','mmflxd');
            dnm = np.copy(dnm[0,region-1,:,:])
            dnm[dnm>1e25] = np.nan
            amoc.setdefault('max',[]).append(np.nanmax(dnm[:,ind1-1:ind2-1])*1e-9)
            amoc.setdefault('26N',[]).append(np.nanmax(dnm[:,ind3])*1e-9)
            print sdate


    time = np.linspace(fyear,lyear,lyear-fyear+1)
    plt.plot(time,amoc['max'],'b',label='max')
    plt.plot(time,amoc['26N'],'r',label='26N')
    plt.legend(loc='lower right')
    return amoc


def amocmean(root_folder,cmpnt,mdl,ext):
    ''' compute mean amoc'''
    prefix,sdate = get_sdate_ini(root_folder,cmpnt,mdl,ext) 
    # set initial values to zero
    amoc, nx , _ ,_ = set_ini_val_zero(prefix, sdate, dim1='region', dim2='depth',
                            dim3='lat') 
    amoc = np.transpose(amoc, (2, 1, 0));
    amoc = timemean(prefix, amoc,'mmflxd',fyear,lyear)
    return amoc


def var3Dmean(root_folder, cmpnt, mdl, ext, varname):
    ''' compute any 3D variable mean'''
    var, nx, ny, nz, = set_ini_val_zero(prefix, sdate, dim1='x', dim2='y',
                                        dim3='depth')
    var = timemean(prefix, var, varname, fyear, lyear)
    return var


def compute_heat_transport(root_folder, project_name, cmpnt, mdl, ext):
    ''' computes ocean and atmosphere heat transport '''

    mw = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float)
    mw = mw/sum(mw)

    prect, nx, ny, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    lHF, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    OLR, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    ASR, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    SHF, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    LWsfc, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    SWsfc, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    Evap, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')
    SnowFlux, _, _, _ = set_ini_val_zero(prefix, sdate, dim1='lon', dim2='lat')

    for year in xrange(np.int(fyear), np.int(lyear)+1):
        for month in xrange(1, 13):
            # control
            filename = root_folder+project_name+'/'+cmpnt+'/hist/'+project_name + \
                     '.'+mdl+'.'+ext+'.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
            #filename = root_folder+project[0]+'/atm/hist/'+project[0] + \
            #         '.cam2.h0.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
            OLR = OLR + ncread_time_surface(filename, 'FLNT', 0, 1, nx, ny) \
                * mw[month-1]
            ASR = ASR + ncread_time_surface(filename, 'FSNT', 0, 1, nx, ny) \
                * mw[month-1]
            LHF = LHF + ncread_time_surface(filename, 'LHFLX', 0, 1, nx, ny) \
                * mw[month-1]
            SHF = SHF + ncread_time_surface(filename, 'SHFLX', 0, 1, nx, ny) \
                * mw[month-1]
            LWsfc = LWsfc + ncread_time_surface(filename, 'FLNS', 0, 1, nx, ny) \
                * mw[month-1]
            SWsfc = SWsfc - ncread_time_surface(filename, 'FSNS', 0, 1, nx, ny) \
                * mw[month-1]
            Evap = Evap + ncread_time_surface(filename, 'QFLX', 0, 1, nx, ny) \
                * mw[month-1]
            precc = ncread_time_surface(filename, 'PRECC', 0, 1, nx, ny) * \
                    mw[month-1]
            precl = ncread_time_surface(filename, 'PRECL', 0, 1, nx, ny) * \
                    mw[month-1]
            precsc = ncread_time_surface(filename, 'PRECSC', 0, 1, nx, ny) * \
                    mw[month-1]
            precsl = ncread_time_surface(filename, 'PRECSL', 0, 1, nx, ny) * \
                    mw[month-1]
            # total precipitation
            prect = prect + (precc+precc)*rhow
            SnowFlux = SnowFlux + (precsc+precsl)*rhow*Lhfus


        print year 

    prect = prect/(np.float(lyear)-np.float(fyear)+1)
    OLR = OLR/(np.float(lyear)-np.float(fyear)+1)
    ASR = ASR/(np.float(lyear)-np.float(fyear)+1)
    LHF = LHF/(np.float(lyear)-np.float(fyear)+1)
    SHF = SHF/(np.float(lyear)-np.float(fyear)+1)
    SnowFlux = SnowFlux/(np.float(lyear)-np.float(fyear)+1)
    LWsfc = LWsfc/(np.float(lyear)-np.float(fyear)+1)
    SWsfc = SWsfc/(np.float(lyear)-np.float(fyear)+1)
    Rtoa = ASR - OLR  # net downwelling radiation
    EminusP = Evap - prect  # kg/m2/s or mm/s
    SurfaceRadiation = LWsfc + SWsfc  # net upward radiation from surface
    SurfaceHeatFlux = SurfaceRadiation + LHF + SHF + SnowFlux  # net upward surface heat flux
    Fatmin = Rtoa + SurfaceHeatFlux  # net heat flux in to atmosphere

    #lon = nc_read(filename, 'lon')
    lat = nc_read(filename, 'lat')
    # mean averaged on x-axis
    Rtoa = np.mean(Rtoa, axis=1)
    SurfaceHeatFlux = np.mean(SurfaceHeatFlux, axis=1) 
    Fatmin = np.mean(Fatmin, axis=1)
    # heat transport terms 
    HTmonthly = {}
    HTmonthly['total'] = inferred_heat_transport(Rtoa, lat)
    HTmonthly['atm'] = inferred_heat_transport(Fatmin, lat)
    HTmonthly['ocn'] = inferred_heat_transport(-SurfaceHeatFlux, lat)
    fig = plt.figure(figsize=(10,4))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(lat_cesm, HT['total'], 'k-', label='total', linewidth=2)
    ax.plot(lat_cesm, HT['atm'], 'r-', label='atm', linewidth=2)
    ax.plot(lat_cesm, HT['ocn'], 'b-', label='ocean', linewidth=2)
    ax.set_xlim(-90,90)
    ax.set_xticks(ticks)
    ax.legend(loc='upper left')
    ax.grid()
    return HTmonthly,lat


def set_ini_val_zero(prefix,sdate,*args, **kwargs):
    dim1 = kwargs.get('dim1', None)
    dim2 = kwargs.get('dim2', None)
    dim3 = kwargs.get('dim3', None)
    dim4 = kwargs.get('dim4', None)
    nx = {}
    ny = {}
    nz = {}
    nx = ncgetdim(prefix+sdate+'.nc', dim1)
    ny = ncgetdim(prefix+sdate+'.nc', dim2)
    var = np.zeros((ny,nx))
    if not dim3 == None:
        nz = ncgetdim(prefix+sdate+'.nc', dim3)
        var = np.zeros((nz,ny,nx))
    if not dim4 == None:
        nk=ncgetdim(prefix+sdate+'.nc', dim4)
        var = np.zeros((nz,ny,nx,nk))


    return var , nx, ny, nz


# 1 for atlantic_arctic_ocean region
# 2 for indian_pacific_ocean region
# 3 for global_ocean
region = 1
root_folder = '/work/milicak/mnt/norstore/NS2345K/noresm/cases/'

expid = 'NOIIA_T62_tn11_norems2_ctrl_tke'
#expid = 'N1850_f19_tn11_01_default'
fyear = 100; # first year
lyear = 110; # last year
m2y = 1
cmpnt = 'ocn' # ocn, atm 
mdl = 'micom' # micom, cam2, cam
ext = 'hm' # hm, hy, h0
#cmpnt = 'atm' # ocn, atm 
#mdl = 'cam2' # micom, cam2, cam
#ext = 'h0' # hm, hy, h0
prefix,sdate = get_sdate_ini(root_folder, cmpnt, mdl, ext) 

#vol=passagevolumetransporttime(root_folder,5) # 5 for Drake
#amoctime=amoctime(root_folder, cmpnt, mdl, ext)
#amocmean=amocmean(root_folder, cmpnt, mdl, ext)
#HT,lat_cesm=compute_heat_transport(root_folder, expid, cmpnt, mdl, ext)
temp=var3Dmean(root_folder, cmpnt, mdl, ext, 'templvl')
