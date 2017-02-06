#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import pdb
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import NorESM_utils as noresmutils

plt.ion()

datesep = '-'


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


def ncread_lev(fname, variable, zlev):
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp = np.squeeze(ncfile.variables[variable][:,zlev,:,:].copy())

    return tmp


def get_sdate_ini(root_folder,cmpnt,mdl,ext, **kwargs):
    expid = kwargs.get('expid', None)
    m2y = kwargs.get('m2y', None)
    # Get dimensions and time attributes for ocn or atm or others ...
    prefix=root_folder+expid+'/'+cmpnt+'/hist/'+expid+ '.'+ mdl + \
            '.'+ext+'.'
    if m2y==1:
        sdate="%4.4d%c%2.2d" % (1,'-',1)    # assuming fyear 1 exists
    else:
        sdate="%4.4d" % (1)


    return prefix, sdate


def inferred_heat_transport(energy_in, lat_deg):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate
    lat_rad = np.deg2rad(lat_deg)
    return (1E-15 * 2 * np.math.pi * aradius**2 *
            integrate.cumtrapz( np.cos(lat_rad)*energy_in, \
            x=lat_rad, initial=0.))


def timemean(prefix, var, varname, fyear, lyear, **kwargs):
    zlev = kwargs.get('zlev', None)
    expid = kwargs.get('expid', None)
    m2y = kwargs.get('m2y', None)
    ''' timemean averages '''
    months2days=[31,  28,  31,  30,  31,   30,   31,  31,   30, 31,   30, 31];
    yeardays=sum(months2days);
    mw = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float)
    mw = mw/sum(mw)
    n=0.0;
    fyear
    lyear
    for year in xrange(fyear,lyear+1):
        if m2y==1:
            for month in xrange(0,12):
                n=n+mw[month]
                sdate="%4.4d%c%2.2d" % (year,'-',month+1)
                if not zlev == None:
                    dnm=ncread_lev(prefix+sdate+'.nc', varname, zlev);
                else:
                    dnm=nc_read(prefix+sdate+'.nc', varname);


                var=var+np.squeeze(dnm.data)*mw[month]
            print sdate, n


        else:
            sdate = "%4.4d" % (year)
            n += 1.0
            if not zlev == None:
                dnm=ncread_lev(prefix+sdate+'.nc', varname, zlev);
            else:
                dnm=nc_read(prefix+sdate+'.nc', varname);


            var=var+np.squeeze(dnm.data)
            print sdate, n


    var = var/n
    return var


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


def zonalmean_bias(root_folder, cmpnt, mdl, ext, varname, woafname, woavname
                  ,maskfile, maskfilevar):
    ''' compute any 3D variable mean'''
    import scipy.io
    var, nx, ny, nz, = set_ini_val_zero(prefix, sdate, dim1='x', dim2='y',
                                        dim3='depth')
    var = timemean(prefix, var, varname, fyear, lyear)
    varwoa = nc_read(woafname, woavname)
    depth = nc_read(woafname, 'depth')
    lat = nc_read(woafname, 'TLAT')
    mat = scipy.io.loadmat(maskfile)
    mask = np.transpose(np.array(mat[maskfilevar]))
    return var, varwoa, mask, lat, depth


def zonalmean_bias_old(root_folder, cmpnt, mdl, ext, varname, woafname, woavname
                  , maskfile, maskindex):
    ''' compute any 3D variable mean'''
    var, nx, ny, nz, = set_ini_val_zero(prefix, sdate, dim1='x', dim2='y',
                                        dim3='depth')
    var = timemean(prefix, var, varname, fyear, lyear)
    mask1 = nc_read(prefix+sdate+'.nc','templvl');
    mask = np.squeeze(mask1.mask)
    print 'mehhmet', mask.shape, mask[0,200,200]
    var_w=ma.zeros((nz,180,360))
    m_w=ma.zeros((nz,180,360))
    #LOOP OVER THE DEPTH LEVELS
    for z in range(var.shape[0]):
        print 'level = ',z
        lon_w,lat_w,t1_w=noresmutils.noresm2WOA(var[z,:-1,:],shift=True,
                                               grid='tnx1v1')
        lon_w,lat_w,m1_w=noresmutils.noresm2WOA(ma.masked_array(1-mask[z,:-1,:],
                            mask=mask[z,:-1,:]),shift=True, grid='tnx1v1')
    # NOTE THAT HERE WE INTERPOLATE THE MASK AS WELL - THIS WILL GIVE US THE
    # LAND FRACTION OF EACH CELL IN THE NEW GRID
    # DIVIDING BY THIS FRACTION GIVES THE TEMPERATURE/SALINITY
    # VALUE OF THE CELL (I.E. NOT THE MEAN OF THE WHOLE CELL,
    # BUT THE MEAN OF THE OCEAN PART)
        var_w[z,:,:]=t1_w
        var_w[z,:,:]=t1_w/m1_w
        m_w[z,:,:]=m1_w


    dnm = np.copy(var_w[:,:,180:])
    var_w[:,:,180:] = var_w[:,:,:180]
    var_w[:,:,:180] = dnm
    return var_w, t1_w, m1_w, mask

def levelvar_bias(root_folder, cmpnt, mdl, ext, varname,woafname,woavname,z1,z2
                 , gridtype):
    ''' compute any 3D variable mean'''
    var, nx, ny, _, = set_ini_val_zero(prefix, sdate, dim1='x', dim2='y',)
    #var, nx, ny, nz, = set_ini_val_zero(prefix, sdate, dim1='x', dim2='y',
    #                                    dim3='depth')
    var = timemean(prefix, var, varname, fyear, lyear, zlev=z1)
    lon_w,lat_w,t1_w=noresmutils.noresm2WOA(var[:-1,:],shift=True,
                                            grid=gridtype)
    # adjust lon_w
    lon_w[lon_w > 180] = lon_w[lon_w > 180]-360
    varsurface = np.copy(t1_w)
    dnm = np.copy(varsurface[:,180:])
    varsurface[:,180:] = varsurface[:,:180]
    varsurface[:,:180] = dnm
    varwoa = nc_read(woafname, woavname)
    varwoa = np.copy(varwoa[z2,:,:])
    lonwoa = nc_read(woafname, 'lon')
    latwoa = nc_read(woafname, 'lat')
    varwoa[varwoa < -10] =0
    return varsurface, varwoa, lonwoa, latwoa


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


def main():
    # tripolar 1degree grid
    grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
    # tripolar 0.25degree grid
    #grid_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/maps/grid_0_25degree.nc';
    # bi-polar grid
    #grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';

    # 1 for atlantic_arctic_ocean region
    # 2 for indian_pacific_ocean region
    # 3 for global_ocean
    region = 1
    root_folder = '/work/milicak/mnt/norstore/NS2345K/noresm/cases/'
    #root_folder = '/work/milicak/mnt/viljework/archive/'

    #expid = 'NOIIA_T62_tn11_norems2_ctrl_tke'
    #expid = 'NBF1850_f19_tn11_sst_sss_rlx_01'
    #expid = 'N1850_f19_tn11_01_default'
    expid = 'NOIIA_T62_tn025_default_visc01'
    fyear = 20; # first year
    lyear = 30; # last year
    m2y = 1
    cmpnt = 'ocn' # ocn, atm
    mdl = 'micom' # micom, cam2, cam
    ext = 'hm' # hm, hy, h0
    varname = 'templvl'
    #cmpnt = 'atm' # ocn, atm
    #mdl = 'cam2' # micom, cam2, cam
    #ext = 'h0' # hm, hy, h0

    gridtype = 'tnx0.25v1' # tnx1v1 , tnx0.25v1, gx1v6
    prefix,sdate = get_sdate_ini(root_folder, cmpnt, mdl, ext)
    woafnamet = '/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/t00an1.nc'
    woafnames = '/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/s00an1.nc'
    woafname = '/fimm/home/bjerknes/milicak/Analysis/obs/WOA13/Analysis/WOA13_tnx1v1_65layers.nc'
    #mask_woa09_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/general/Analysis/woa_mask.mat';
    mask_woa09_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/general/Analysis/noresm_tnxv1_mask.mat';

    #mask_index = 0; # 0 for Global
    mask_index = 10; # 10 for Atlantic Ocean
    #mask_index = 1; # 1 for Arctic Ocean
    #mask_index = 2; # 2 for Mediterranean
    #mask_index = 3; # 3 for Pacific Ocean
    #mask_index = 4; # 4 for Southern Ocean
    #mask_index = 6; # 6 for Baltic Sea
    #mask_index = 7; # 7 for Red Sea
    #mask_index = 8; # 8 for Indian Ocean
    #mask_index = 9; # 9 for Black Sea and Caspian Sea


    def call_generic_diags(argument):
        switcher = {
            'amoctime': 1,
            'amocmean': 2,
            'heattransport': 3,
            'voltr': 4,
            '3Dmean': 5,
            'sstbias': 6,
            'sssbias': 7,
            'zonalmean': 8,
        }
        return switcher.get(argument, "non-valid option. Please select another option")


    diagno = call_generic_diags('sstbias')
    print 'You select', diagno

    if diagno == 1:
        amoctime = amoctime(root_folder, cmpnt, mdl, ext)
    elif diagno == 2:
        amocmean = amocmean(root_folder, cmpnt, mdl, ext)
    elif diagno == 3:
        HT,lat_cesm = compute_heat_transport(root_folder, expid, cmpnt, mdl, ext)
    elif diagno == 4:
        vol = passagevolumetransporttime(root_folder,5) # 5 for Drake
    elif diagno == 5:
        temp = var3Dmean(root_folder, cmpnt, mdl, ext, varname)
    elif diagno == 6:
        sst,sstwoa,lon,lat = levelvar_bias(root_folder, cmpnt, mdl, ext, 'templvl',
                                      woafnamet,'t',0,0,gridtype)
        plt.figure()
        m=Basemap(llcrnrlon=0,llcrnrlat=-88,urcrnrlon=360,urcrnrlat=88,projection='cyl')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
        m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
        im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sst-sstwoa),
                           shading='flat',vmin=-5,vmax=5,cmap='RdBu_r');
        cb = m.colorbar(im1,"right", size="5%", pad="10%")
        #plt.pcolor(lon,lat,np.ma.masked_invalid(sst-sstwoa),vmin=-5,vmax=5);plt.colorbar()
    elif diagno == 7:
        sss,ssswoa,lon,lat = levelvar_bias(root_folder, cmpnt, mdl, ext, 'salnlvl',
                                      woafnames,'s',0,0,gridtype)
        plt.figure()
        m=Basemap(llcrnrlon=0,llcrnrlat=-88,urcrnrlon=360,urcrnrlat=88,projection='cyl')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
        m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
        im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sss-ssswoa),
                           shading='flat',vmin=-3,vmax=3,cmap='RdBu_r');
        cb = m.colorbar(im1,"right", size="5%", pad="10%")
    elif diagno == 8:
        var, varwoa, mask, lat, depth = zonalmean_bias(root_folder, cmpnt, mdl, ext,
                                                       'templvl',
                                      woafname, 'twoa_noresm', mask_woa09_file,
                                                 'mask_tnxv1')
        tmask = np.copy(var)
        tmask[tmask>-50.0]=1.0
        tmask[tmask<=-50.0]=0.0
        var[var<=-50]=0.
        mask = np.double(mask)
        if mask_index == 0:
            mask[mask != mask_index] = 1.0
            mask[mask == mask_index] = np.nan
        else:
            mask[mask != mask_index] = np.nan
            mask[mask == mask_index] = 1.0


        area = nc_read(grid_file, 'parea')
        zonalbias = np.zeros((varwoa.shape[0],area.shape[0],area.shape[1]))
        zonalbiaswght = np.zeros((varwoa.shape[0],area.shape[0],area.shape[1]))
        for z in xrange(0,varwoa.shape[0]):
            zonalbias[z,:,:] = np.squeeze(var[z,:,:]-varwoa[z,:,:])*mask*area
            zonalbiaswght[z,:,:] = np.squeeze(tmask[z,:,:])*mask*area


        zonalbias = np.nansum(zonalbias, axis=2)
        zonalbiaswght = np.nansum(zonalbiaswght, axis=2)
        zonalbias = zonalbias/zonalbiaswght
        plt.figure()
        plt.pcolor(lat[:,0],-depth,np.ma.masked_invalid(zonalbias),vmin=-5,vmax=5);plt.colorbar()
        #plt.pcolor(lat,depth,np.ma.masked_invalid(var-varwoa),vmin=-5,vmax=5);plt.colorbar()




if __name__ == "__main__":
    # this won't be run when imported
    main()

