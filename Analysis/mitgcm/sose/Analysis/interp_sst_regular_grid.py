import numpy as np
import pandas as pd
from vcr import utils, conserve
from HCtFlood.kara import flood_kara
import xesmf as xe

def convert_lon(ds, lon='XC'):
    ds.coords[lon] = (ds.coords[lon] + 180) % 360 - 180
    ds = ds.sortby(ds[lon])
    return ds


root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'

grname = root_folder+project_name+'/grid.nc'
woa = xr.open_dataset('oisst-avhrr-v02r01.20060801.nc',decode_times=False)
mask = xr.where(woa.sst[0,0,:,:]>-10,1,0)
oi = woa
oi = oi.drop({'sst', 'anom', 'err', 'ice', 'time', 'zlev'})
gr = xr.open_dataset(grname)
gr = gr.drop_dims('time')
# encoding={'sst': {'dtype': 'float32', 'scale_factor': 0.01, '_FillValue':
#                   -999},
#           'ice': {'dtype': 'float32', 'scale_factor': 0.01, '_FillValue':
#                   -999}}
encoding={'sst': {'dtype': 'int16', 'scale_factor': 0.01, '_FillValue':
                  -999},
          'ice': {'dtype': 'int16', 'scale_factor': 0.01, '_FillValue':
                  -999}}


# for year in range(2005,2013):
for year in range(2011,2013):
    for month in range(1,13):
        date = np.str(year) + '-' + np.str(month).zfill(2)
        p = pd.Period(date)
        days = p.days_in_month
        for ddays in range(1,days+1):
            print(date+'-'+np.str(ddays))
            fnames = root_folder+project_name+'/SST_' + np.str(year) + '_' + np.str(month).zfill(2) + '_' + np.str(ddays).zfill(2) +'.nc'
            df = xr.open_dataset(fnames)
            fnames = root_folder+project_name+'/SIarea_' + np.str(year) + '_' + np.str(month).zfill(2) + '_' + np.str(ddays).zfill(2) +'.nc'
            df2 = xr.open_dataset(fnames)
            df = xr.merge([df,df2])
            df = df.rename_dims({'i': 'XC', 'i_g': 'XG', 'j': 'YC', 'j_g': 'YG'})
            df = df.rename_dims({'k': 'Z', 'k_u': 'ZU', 'k_l': 'Zl', 'k_p1': 'Zp1'})
            ds = xr.merge([df, gr])
            ds = ds.drop({'maskC', 'maskW', 'maskS', 'maskInC', 'maskInW', 'maskInS',
                          'hFacC', 'hFacW', 'hFacS', 'dxG', 'dyG', 'dxC', 'dyC', 'rAw',
                          'rAs'})
            # convert lon between -180 to 180
            # ds = convert_lon(ds)

            ds1 = ds
            ds1 = ds1.drop({'i', 'i_g', 'j', 'j_g', 'k', 'k_u', 'k_l',
                    'rA', 'Depth', 'THETA', 'rAz', 'XG', 'YG', 'Z', 'Zl', 'Zu'})

            regrid_domain = xe.Regridder(ds1.rename({'XC': 'lon', 'YC': 'lat'}),
                             oi, 'nearest_s2d',
                             periodic=True,
                             filename='/archive2/milicak/mitgcm/sose/regrid_domain.nc',
                             reuse_weights=True)


            # extrapolate ocean values into the land
            # drowned_temp = flood_kara(ds['THETA'], xdim='XC', ydim='YC')
            # interpoalte to the regional domain
            sst_woa = regrid_domain(ds.THETA)
            ice_woa = regrid_domain(df2.SIarea)
            # replace 0 with nan
            mask = xr.where(woa.sst[0,0,:,:]>-10,1,0)
            mask = mask.where(mask!=0)
            sst_woa = sst_woa*mask
            ice_woa = ice_woa*mask
            # ice_woa = ice_woa*100
            sst = sst_woa.to_dataset(name='sst')
            ice = ice_woa.to_dataset(name='ice')
            dsnew = xr.merge([sst, ice])
            dsnew.lon.attrs['units'] = 'degrees_east'
            dsnew.lat.attrs['units'] = 'degrees_north'
            dsnew.lon.attrs['grids'] = woa.lon.grids
            dsnew.lat.attrs['grids'] = woa.lat.grids
            outputfile = root_folder+project_name+'/WOA_grid_SST_SI_' + np.str(year) + '_' + np.str(month).zfill(2) + '_' + np.str(ddays).zfill(2) +'.nc'
            dsnew.to_netcdf(outputfile,encoding=encoding)





# outputfile1 = root_folder+project_name+'/SST_2011_01_08.bin'
# outputfile2 = root_folder+project_name+'/SIarea_2011_01_08.bin'
# outputfile3 = root_folder+project_name+'/SST_header_2011_01_08.bin'
# sst_woa.data.astype('>f4').tofile(outputfile1)
# ice_woa.data.astype('>f4').tofile(outputfile2)
# aa = np.array([2011,1,8])
# aa.astype('>f4').tofile(outputfile3)



# polar projection
# delx = 1.0/12
# dely = 1.0/60
# lon = np.arange(0.5*delx,360, delx)
# lat = np.arange(-78.0,0, dely)
#
# dsnew = xr.Dataset(data_vars={"lon":(["lon"],lon),
#                            "lat":(["lat"],lat)})
#
# regrid_domain = xe.Regridder(ds1.rename({'XC': 'lon', 'YC': 'lat'}),
#                              dsnew
#                              , 'nearest_s2d',
#                             periodic=True,
#                             filename='/archive2/milicak/mitgcm/sose/regrid_domain_polar.nc',
#                             reuse_weights=True)
#
#
