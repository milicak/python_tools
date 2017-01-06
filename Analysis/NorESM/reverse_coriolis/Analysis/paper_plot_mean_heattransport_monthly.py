'''
This routine computes the ocean and atm heat transport annually. Inspired from:
http://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes
/Lectures/Lecture13%20--%20Heat%20transport.html#section6
'''
import my_nanfilter
import numpy as np
import numpy.ma as ma
# import scipy.io
# from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from netCDF4 import Dataset
from netcdf_functions import nc_read
import matplotlib.pyplot as plt
import matplotlib as mpllib
import sys
reload(my_nanfilter)

# from my_nanfilter import my_nanfilterbox
# from netcdf_functions import ncgetdim
# %matplotlib inline
# np.shape !!!!!

# IMPORTANT 
plt.ion()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)


def inferred_heat_transport(energy_in, lat_deg):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate
    lat_rad = np.deg2rad(lat_deg)
    aradius = 6.373E6      # Radius of Earth (m)
    return (1E-15 * 2 * np.math.pi * aradius**2 *
            integrate.cumtrapz( np.cos(lat_rad)*energy_in, \
            x=lat_rad, initial=0.))

def ncread_time_surface(fname, variable, timestr, timeend, x, y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp = np.zeros([y, x])
    for i in range(timestr, timeend):
        # print i
        tmp = tmp+ncfile.variables[variable][i, :, :].copy()

    tmp = tmp/(timeend-timestr)
    return tmp


def enable_global(tlon,tlat,data):
    """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
    tlon = np.where(np.greater_equal(tlon, min(tlon[:, 0])), tlon-360, tlon)
    tlon = tlon+abs(ma.max(tlon))
    tlon = tlon+360
    # stack grids side-by-side (in longitiudinal direction), so
    # any range of longitudes may be plotted on a world map.
    tlon = np.concatenate((tlon, tlon+360), 1)
    tlat = np.concatenate((tlat, tlat), 1)
    data = ma.concatenate((data, data), 1)
    tlon = tlon-360.
    return tlon, tlat, data



def compute_heat_transport(root_folder,project_name,cam_ext):
    ''' computes heat transport '''
    Lv = 2.5e6  # Joule/kg
    rhow = 1000  # kg/m3
    Lhvap = 2.5E6    # Latent heat of vaporization (J / kg)
    Lhsub = 2.834E6   # Latent heat of sublimation (J / kg)
    Lhfus = Lhsub - Lhvap  # Latent heat of fusion (J / kg)

    mw = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float)
    mw = mw/sum(mw)

    prect = np.zeros([ny, nx])
    LHF = np.zeros([ny, nx])
    OLR = np.zeros([ny, nx])
    ASR = np.zeros([ny, nx])
    SHF = np.zeros([ny, nx])
    LWsfc = np.zeros([ny, nx])
    SWsfc = np.zeros([ny, nx])
    Evap = np.zeros([ny, nx])
    SnowFlux = np.zeros([ny, nx])

    #sys.exit()
    for year in xrange(np.int(fyear), np.int(lyear)+1):
        for month in xrange(1, 13):
            # control
            filename = root_folder+project_name[0]+'/atm/hist/'+project_name[0] + \
                     '.'+cam_ext+'.h0.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
            #filename = root_folder+project[0]+'/atm/hist/'+project[0] + \
            #         '.cam2.h0.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
            # reverse Coriolis
            #filename = root_folder+project_rf[0]+'/atm/hist/'+project_rf[0] + \
            #        '.cam.h0.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
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
    return HTmonthly,lat


root_folder = '/work/milicak/mnt/norstore/NS2345K/noresm/cases/'
project = ['N1850_f19_tn11_01_default']
project_rf = ['N1850_f19_tn11_reverseCoriolis']
titles = ['ctrl','reverse f']

fyear = '650'  # first year
lyear = '750'  # last year
nx = 144
ny = 96
ticks = [-90, -60, -30, 0, 30, 60, 90]


fig = plt.figure(figsize=(10,4))
runs = [project, project_rf]
cams = ['cam2','cam']
N = len(runs)

for n, HTname in enumerate([project, project_rf]):
    HT,lat_cesm=compute_heat_transport(root_folder,HTname,cams[n])
    ax = fig.add_subplot(1, N, n+1)
    ax.plot(lat_cesm, HT['total'], 'k-', label='total', linewidth=2)
    ax.plot(lat_cesm, HT['atm'], 'r-', label='atm', linewidth=2)
    ax.plot(lat_cesm, HT['ocn'], 'b-', label='ocean', linewidth=2)
    ax.set_title(titles[n])
    ax.set_xlim(-90,90)
    ax.set_xticks(ticks)
    ax.legend(loc='upper left')
    ax.grid()


plt.savefig('paperfigs/total_heat_transport.eps',
            bbox_inches='tight', format='eps', dpi=200)
plt.clf()
plt.close(fig)
sys.exit()

fig = plt.figure()
# ax = plt.gca()
# ax.set_axis_bgcolor('grey')
# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
# m = Basemap(projection='kav7',lon_0=0,resolution='c')
# m = Basemap(llcrnrlon=-180, llcrnrlat=-88, urcrnrlon=180, urcrnrlat=88,
#             projection='cyl')
# m.drawcoastlines()
# m.drawparallels(np.arange(-80, 81, 20), labels=[1, 1, 0, 0])
# m.drawmeridians(np.arange(0, 360, 60), labels=[0, 0, 0, 1])
# im1 = m.pcolormesh(lon, lat, np.ma.masked_invalid(EminusP), shading='flat',
#                    cmap='jet', vmin=-13, vmax=5, latlon=True)
# cb = m.colorbar(im1, "right", size="5%", pad="10%")
# cb.set_label('[mm/day]')
# plt.ylabel('Depth [m]')
# plt.xlabel('Lat')
# plt.ylim(-7000,0)
# cb.set_label('[' r'$^\circ$' 'C]')
plt.savefig('paperfigs/'+project[0]+'_heat_transport.eps',
            bbox_inches='tight', format='eps', dpi=200)
plt.clf()
plt.close(fig)


# plt.show()



