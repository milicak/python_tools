#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
import pdb
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
#import NorESM_utils as noresmutils
import argparse
import glob
import ESMF

plt.ion()

datesep = '-'


# global constants
aradius = 6.373E6      # Radius of Earth (m)
Lv = 2.5e6  # Joule/kg
rhow = 1000  # kg/m3
Lhvap = 2.5E6    # Latent heat of vaporization (J / kg)
Lhsub = 2.834E6   # Latent heat of sublimation (J / kg)
Lhfus = Lhsub - Lhvap  # Latent heat of fusion (J / kg)



def get_sdate_ini(root_folder,cmpnt,mdl,ext, **kwargs):
    expid = kwargs.get('expid', None)
    m2y = kwargs.get('m2y', None)
    # Get dimensions and time attributes for ocn or atm or others ...
    #prefix=root_folder+expid+'/'+'/run/'+expid+ '.'+ mdl + \
    #        '.'+ext+'.'
    prefix=root_folder+expid+'/'+cmpnt+'/hist/'+expid+ '.'+ mdl + \
            '.'+ext+'.'
    if m2y==1:
        sdate="%4.4d%c%2.2d" % (1,'-',1)    # assuming fyear 1 exists
    else:
        sdate="%4.4d" % (1)


    return prefix, sdate

def amocmean(root_folder, expid, cmpnt, mdl, ext, region):
    ''' compute max amoc, amoc26 in time'''

    # Latitude interval for maximum AMOC search
    foldername = root_folder + '/' + expid + '/' + cmpnt + '/hist/'
    sdate="%c%4.4d%c" % ('*',fyear,'*')
    freq = '*'+ext
    list = sorted(glob.glob(foldername+freq+sdate))
    for year in xrange(fyear+1,lyear+1):
        sdate = "%c%4.4d%c" % ('*',year,'*')
        list.extend(sorted(glob.glob(foldername+freq+sdate)))

    chunks = (385,360)
    xr_chunks = {'x': chunks[-2], 'y': chunks[-2]}
    data = xr.open_mfdataset(list,decode_times=False, chunks=xr_chunks)['mmflxd']
    dnm = np.copy(data[:,region-1,:,:].mean('time'))*1e-9
    plt.figure()
    plt.pcolor(data.lat,-data.depth,np.ma.masked_invalid(dnm),
               vmin=-10,vmax=30);plt.colorbar()
    return dnm

def amoctime(root_folder, expid, cmpnt, mdl, ext, region):
    ''' compute max amoc, amoc26 in time'''

    # Latitude interval for maximum AMOC search
    foldername = root_folder + '/' + expid + '/' + cmpnt + '/hist/'
    sdate="%c%4.4d%c" % ('*',fyear,'*')
    freq = '*'+ext
    list = sorted(glob.glob(foldername+freq+sdate))

    chunks = (385,360)
    xr_chunks = {'x': chunks[-2], 'y': chunks[-2]}
    data = xr.open_mfdataset(list,decode_times=False, chunks=xr_chunks)['mmflxd']
    lat1 = 20;
    lat2 = 60;
    lat3 = 26.5;
    ind1 = np.min(np.where(data.lat>=lat1));
    ind2 = np.max(np.where(data.lat<=lat2));
    ind3 = np.max(np.where(data.lat<=lat3));
    amoc = {}
    dnm = np.nanmax(np.copy(data[:,region-1,:,ind1-1:ind2-1].mean('time')))*1e-9
    amoc.setdefault('max',[]).append(dnm)
    dnm = np.nanmax(np.copy(data[:,region-1,:,ind3-1].mean('time')))*1e-9
    amoc.setdefault('26N',[]).append(dnm)
    for year in xrange(fyear+1,lyear+1):
        sdate = "%c%4.4d%c" % ('*',year,'*')
        list = sorted(glob.glob(foldername+freq+sdate))
        data = xr.open_mfdataset(list,decode_times=False, chunks=xr_chunks)['mmflxd']
        dnm = np.nanmax(np.copy(data[:,region-1,:,ind1-1:ind2-1].mean('time')))*1e-9
        amoc.setdefault('max',[]).append(dnm)
        dnm = np.nanmax(np.copy(data[:,region-1,:,ind3-1].mean('time')))*1e-9
        amoc.setdefault('26N',[]).append(dnm)


    time = np.linspace(fyear,lyear,lyear-fyear+1)
    plt.figure()
    line1, = plt.plot(time,amoc['max'],'k',label='max')
    line2, = plt.plot(time,amoc['26N'],'r',label='26N')
    plt.legend(loc='lower right')
    return amoc, line1, line2

def passagevolumetransporttime(root_folder, expid, cmpnt, mdl, ext, region):
    ''' compute Drake Passage transport in time'''

    # Latitude interval for maximum AMOC search
    foldername = root_folder + '/' + expid + '/' + cmpnt + '/hist/'
    sdate="%c%4.4d%c" % ('*',fyear,'*')
    freq = '*'+ext
    list = sorted(glob.glob(foldername+freq+sdate))

    chunks = (385,360)
    xr_chunks = {'x': chunks[-2], 'y': chunks[-2]}
    data = xr.open_mfdataset(list,decode_times=False,
                             chunks=xr_chunks)['voltr']
    vol= {}
    dnm = np.copy(data[:,region-1].mean('time'))*1e-9
    vol.setdefault('Drake',[]).append(dnm)
    for year in xrange(fyear+1,lyear+1):
        sdate = "%c%4.4d%c" % ('*',year,'*')
        list = sorted(glob.glob(foldername+freq+sdate))
        data = xr.open_mfdataset(list,decode_times=False,
                                 chunks=xr_chunks)['voltr']
        dnm = np.copy(data[:,region-1].mean('time'))*1e-9
        vol.setdefault('Drake',[]).append(dnm)


    time = np.linspace(fyear,lyear,lyear-fyear+1)
    plt.figure()
    line1, = plt.plot(time,vol['Drake'],'k',label='Drake')
    plt.legend(loc='lower right')
    return vol


def enable_global(tlon,tlat,data):
  """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
  tlon = tlon+abs(ma.max(tlon)); tlon=tlon+360
  # stack grids side-by-side (in longitiudinal direction), so
  # any range of longitudes may be plotted on a world map.
  tlon = np.concatenate((tlon,tlon+360),1)
  tlat = np.concatenate((tlat,tlat),1)
  data = ma.concatenate((data,data),1)
  tlon = tlon-360.
  return tlon, tlat, data


def levelvar_bias_xr(root_folder, cmpnt, mdl, ext, varname,woafname,woavname,z1,z2
                 , gridtype, **kwargs):
    ''' compute any 3D variable mean'''
    prefix = kwargs.get('prefix', None)
    sdate = kwargs.get('sdate', None)
    m2y = kwargs.get('m2y', None)
    expid = kwargs.get('expid', None)
    foldername = root_folder + '/' + expid + '/' + cmpnt + '/hist/'
    sdate="%c%4.4d%c" % ('*',fyear,'*')
    freq = '*'+ext
    list = sorted(glob.glob(foldername+freq+sdate))
    for year in xrange(fyear+1,lyear+1):
        sdate = "%c%4.4d%c" % ('*',year,'*')
        list.extend(sorted(glob.glob(foldername+freq+sdate)))


    chunks = (385,360)
    xr_chunks = {'x': chunks[-2], 'y': chunks[-2]}
    data = xr.open_mfdataset(list,decode_times=False, chunks=xr_chunks)[varname]
    var = np.copy(data[:,z1,:,:].mean('time'))
    # create source grid from noresm grid to WOA grid/data netcdf file
    srcgrid = ESMF.Grid(filename=grid_file, filetype=ESMF.FileFormat.GRIDSPEC,
                       coord_names=["vlon","vlat"])
    # create a destination grid file
    dstgrid = ESMF.Grid(filename=woafname, filetype=ESMF.FileFormat.GRIDSPEC)
    # Create a field on the centers of the grid
    field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
    field1.data[:] = np.transpose(var)
    # Create a field on the centers of the grid
    field2 = ESMF.Field(dstgrid, staggerloc=ESMF.StaggerLoc.CENTER)
    # set it to zero
    field2.data[:] = 0
    # Set up a regridding object between source and destination
    regridS2D = ESMF.Regrid(field1, field2,
                            regrid_method=ESMF.RegridMethod.BILINEAR)
    field2 = regridS2D(field1, field2)
    varsurface = field2.data

    varwoa = xr.open_dataset(woafname,decode_times=False)[woavname]
    lon = np.copy(varwoa.lon)
    lat = np.copy(varwoa.lat)
    varwoa = np.copy(varwoa[0,z2,:,:])
    varwoa = np.transpose(varwoa)
    varwoa[varwoa < -10] =0
    #return varsurface, varwoa, lonwoa, latwoa
    return varsurface, varwoa, lon, lat

def main():
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


    global fyear, lyear
    global grid_file
    global lon, lat
    print 'input order = root_folder expid fyear lyear m2y cmpnt mdl ext varname diagname'
    # general_diagnostics.py /work/milicak/mnt/norstore/NS2345K/noresm/cases/
    # NOIIA_T62_tn11_FAMOS_BG_CTR 11 15 ocn micom hm templvl tnx1v1 sstbias

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

    # 1 for atlantic_arctic_ocean region
    # 2 for indian_pacific_ocean region
    # 3 for global_ocean
    region = 1

    root_folder = str(sys.argv[1]) # root_folder
    expid = str(sys.argv[2]) # which experiment
    fyear = int(sys.argv[3]) # first year
    lyear = int(sys.argv[4]) #last year
    cmpnt = str(sys.argv[5])  # component = ocn, atm, cice
    mdl = str(sys.argv[6]) # model micom, pop, cam2, camp
    ext = str(sys.argv[7]) # extention, hm, hy, h0
    varname = str(sys.argv[8]) # variable name sst, templvl, salnlvl
    gridtype = str(sys.argv[9]) #grid_type tnx1v1, tnxv025v1
    diagname = str(sys.argv[10]) #diagnostic name
    print 'expid = ', expid

    #gridtype = 'tnx1v1' # tnx1v1 , tnx0.25v1, gx1v6
    print gridtype
    if gridtype=='tnx1v1':
        # tripolar 1degree grid
        mask_woa09_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/general/Analysis/noresm_tnxv1_mask.mat';
        maskvariable = 'mask_tnxv1'
        grid_file = '/tos-project1/NS2345K/noresm/inputdata/ocn/micom/tnx1v1/20120120/grid.nc';
        woafnamet = 'tempwoa_noresm1deg.nc'
        woafnames = 'saltwoa_noresm1deg.nc'
    elif gridtype=='tnx0.25v1':
        # tripolar 0.25degree grid
        mask_woa09_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/general/Analysis/noresm_tnx0_25v1_mask.mat';
        maskvariable = 'mask'
        grid_file = '/tos-project1/NS2345K/noresm/inputdata/ocn/micom/tnx0.25v4/20170619/grid.nc'
        #woafnamet = 'tempwoa_noresm0_25deg.nc'
        woafnamet = '/tos-project1/NS2345K/noresm_diagnostics/packages/MICOM_DIAG/obs_data/WOA13/0.25deg/woa13_decav_t00_04.nc'
        woafnames = '/tos-project1/NS2345K/noresm_diagnostics/packages/MICOM_DIAG/obs_data/WOA13/0.25deg/woa13_decav_s00_04.nc'
    # bi-polar grid
    #grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';




    #diagname = 'sstbias' #ssstbias, sssbias, zonalmean
    print 'You select', diagname
    diagno = call_generic_diags(diagname)

    if diagno == 1:
        global amoc_time,line1, line2
        amoc_time,line1,line2 = amoctime(root_folder, expid, cmpnt, mdl, ext, region)
        # to change label names later
        # plt.legend([line1, line2], ['Line Up', 'Line Down'])
    elif diagno == 2:
        global amoc_mean
        amoc_mean = amocmean(root_folder, expid, cmpnt, mdl, ext, region)
    elif diagno == 3:
        HT,lat_cesm = compute_heat_transport(root_folder, expid, cmpnt, mdl, ext)
    elif diagno == 4:
        vol = passagevolumetransporttime(root_folder, expid, cmpnt, mdl, ext,
                                         5) # 5 for Drake
    elif diagno == 5:
        temp = var3Dmean(root_folder, cmpnt, mdl, ext, varname)
    elif diagno == 6:
        global sst, sstwoa
        sst,sstwoa,lon,lat = levelvar_bias_xr(root_folder, cmpnt, mdl, ext, 'templvl',
                                      woafnamet,'t_an',0,0,gridtype,
                                             expid=expid)
        plt.figure()
        #m = Basemap(llcrnrlon=280,llcrnrlat=20,urcrnrlon=360,urcrnrlat=80,projection='cyl')
        m = Basemap(llcrnrlon=-180,llcrnrlat=-88,urcrnrlon=180,urcrnrlat=90,projection='cyl')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
        m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
        #m.drawparallels(np.arange(20,80,10),labels=[1,1,0,0])
        #m.drawmeridians(np.arange(280,360,10),labels=[0,0,0,1])
        im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(sst-sstwoa)),
                           shading='flat',vmin=-3,vmax=3,cmap='RdBu_r');
        cb = m.colorbar(im1,"right", size="5%", pad="10%")
        #plt.pcolor(lon,lat,np.ma.masked_invalid(sst-sstwoa),vmin=-5,vmax=5);plt.colorbar()
    elif diagno == 7:
        global sss, ssswoa
        sss,ssswoa,lon,lat = levelvar_bias_xr(root_folder, cmpnt, mdl, ext, 'salnlvl',
                                      woafnames,'s_an',0,0,gridtype,
                                          expid=expid)
        plt.figure()
        m = Basemap(llcrnrlon=-180,llcrnrlat=-88,urcrnrlon=180,urcrnrlat=90,projection='cyl')
        #m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
        m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
        im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(sss-ssswoa)),
                           shading='flat',vmin=-3,vmax=3,cmap='RdBu_r');
        cb = m.colorbar(im1,"right", size="5%", pad="10%")
    elif diagno == 8:
        global zonalbias, depth
        var, varwoa, mask, lat, depth = zonalmean_bias(root_folder, cmpnt, mdl, ext,
                                                       'templvl',
                                      woafname, 'twoa_noresm', mask_woa09_file,
                                                 maskvariable,
                                          prefix=prefix, sdate=sdate, m2y=m2y)
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
            zonalbias[z,:,:] = np.squeeze(var[z,:,:]-varwoa[z,:,:])*area*mask
            zonalbiaswght[z,:,:] = np.squeeze(tmask[z,:,:])*area*mask


        zonalbias = np.nansum(zonalbias, axis=2)
        zonalbiaswght = np.nansum(zonalbiaswght, axis=2)
        zonalbias = zonalbias/zonalbiaswght
        plt.figure()
        plt.pcolor(lat[:,0],-depth,np.ma.masked_invalid(zonalbias),vmin=-3,vmax=3);plt.colorbar()
        plt.set_cmap('RdBu_r')
        plt.xlim(-40,60)
        #plt.pcolor(lat,depth,np.ma.masked_invalid(var-varwoa),vmin=-5,vmax=5);plt.colorbar()


if __name__ == "__main__":
    # this won't be run when imported
    main()

