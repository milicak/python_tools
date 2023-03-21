import numpy as np
import xarray as xr
import glob

#    sn_wndi = 'drowned_uas_JRA55'    ,  3.         ,  'uas'    ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bicubic.nc' , 'U1' ,   ''
#    sn_wndj = 'drowned_vas_JRA55'    ,  3.         ,  'vas'    ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bicubic.nc' , 'V1' ,   ''
#    sn_qsr  = 'drowned_rsds_JRA55'   ,  3.         ,  'rsds'   ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''
#    sn_qlw  = 'drowned_rlds_JRA55'   ,  3.         ,  'rlds'   ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''
#    sn_tair = 'drowned_tas_JRA55'    ,  3.         ,  'tas'    ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''
#    sn_humi = 'drowned_huss_JRA55'   ,  3.         ,  'huss'   ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''
#    sn_prec = 'drowned_tprecip_JRA55' , 3.         ,  'tprecip',    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''
#    sn_snow = 'drowned_prsn_JRA55 '  ,  3.         ,  'prsn'   ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''
#    sn_slp  = 'drowned_psl_JRA55'    ,  3.         ,  'psl'    ,    .true.   , .false. , 'yearly'  , 'eORCA025_JRA55_do_c3.0_weights_bilin.nc' , ''   ,   ''

# Functions for humidity borrowed and adapted from MetPy.calc: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.html
def mixing_ratio(partial_press, total_press, molecular_weight_ratio=0.622):     
    return (molecular_weight_ratio * partial_press                              
                / (total_press - partial_press))                                
                                                                                
                                                                                
def specific_humidity_from_mixing_ratio(mr):                                    
    return mr / (1 + mr)                                                        
                                                                                
                                                                                
def saturation_vapor_pressure(temperature):                                     
    sat_pressure_0c = 6.112e2 # Pa                                              
    return sat_pressure_0c * np.exp(17.67 * (temperature - 273.15) # K -> C     
                                        / (temperature - 29.65))   # K -> C     
                                                                                
def saturation_mixing_ratio(total_press, temperature):                          
    return mixing_ratio(saturation_vapor_pressure(temperature), total_press)    
                                                                                


root_folder = '/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/'
out_folder = '/work/opa/mi19918/Projects/nemo/dataset/ECMWF/'

# -d lon,-19.,42. -d lat,30.,48.

rho0 = 1e3
pprec_fqh = 6
year = 2019
# ls1 = sorted(glob.glob(root_folder + str(year) + '/*/*-MEDATL*fv10*.nc'))
# ls2 = ls1[1::2]
# ls1 = ls1[0::2]
ls1 = sorted(glob.glob(root_folder + str(year) + '/*/*BLK*.nc'))

for idx, fname in enumerate(ls1):
    df = xr.open_dataset(fname)
    df = df.reindex(lat=list(reversed(df.lat))) 
    df = df.rename({'time': 'time_counter','lat': 'y', 'lon': 'x'})
    # convert variables
    # total precip m to kg/m2/s
    df['TP'] = df.TP*rho0/(pprec_fqh*3600)
    df['TP'].attrs['units'] = 'kg/m2/s' 
    # solid precip (snow) m to kg/m2/s
    df['SF'] = df.SF*rho0/(pprec_fqh*3600)
    df['SF'].attrs['units'] = 'kg/m2/s' 
    # msl pa to hPa
    # df['MSL'] = df.MSL*0.01
    # df['MSL'].attrs['units'] = 'hPa' 
    # convert radiation from J/m2 to W/m2: https://confluence.ecmwf.int/pages/viewpage.action?pageId=155337784
    df['SSRD'] = df.SSRD/(pprec_fqh*3600)
    # df['SSRD'] = 80*df.SSRD/df.SSRD
    df['SSRD'].attrs['units'] = 'W m-2' 
    df['STRD'] = df.STRD/(pprec_fqh*3600)
    # df['STRD'] = 320*df.STRD/df.STRD
    df['STRD'].attrs['units'] = 'W m-2' 
    smr = saturation_mixing_ratio(df.SP, df.D2M)          
    sphum = specific_humidity_from_mixing_ratio(smr)   
    sphum.name = 'huss'        
    sphum = sphum.to_dataset() 
    df['huss'] = sphum.huss
    df = df.drop_vars(('SP','D2M'))
    # format for output ecmwf_y2013m03d25.nc
    outname = out_folder + 'ecmwf_y' + str(year) + 'm' + fname[-47:-45] + 'd' + fname[-45:-43]  + '.nc'
    # Remove all _FillValue                                             
    all_vars = list(df.data_vars.keys()) + list(df.coords.keys()) 
    encodings = {v: {'_FillValue': None} for v in all_vars}             
    # Also fix the time encoding                                        
    # encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'hours since 1900-01-01 00:00:00'})
    encodings['time_counter'].update({'dtype':'float64'})
    print(outname)
    df.to_netcdf(                 
    outname,                        
    format='NETCDF4_CLASSIC',    
    engine='netcdf4',            
    encoding=encodings,          
    unlimited_dims=['time_counter'])                                
    df.close()                    



# for 2013 12 31
if year == 2014:
    dd = xr.open_dataset('/data/inputs/metocean/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2013/12/20131231-ECMWF---AM0125-MEDATL-b20140101_an-fv10.00.nc')
    dd = dd.rename({'time': 'time_counter'})
    df = xr.open_dataset('/work/opa/mi19918/Projects/nemo/dataset/ECMWF/ecmwf_y2014m01d01.nc')
    df['time_counter'] = dd.time_counter
    outname = '/work/opa/mi19918/Projects/nemo/dataset/ECMWF/ecmwf_y2013m12d31.nc'
    df.to_netcdf(                 
    outname,                        
    format='NETCDF4_CLASSIC',    
    engine='netcdf4',            
    encoding=encodings,          
    unlimited_dims=['time_counter'])                                
    df.close()                    


# for idx, fname in enumerate(ls1):
#     df = xr.open_dataset(fname)
#     df = df.isel(lon=slice(448,938),lat=slice(176,321))
    # df = df.reindex(lat=list(reversed(df.lat))) 
    # df = df.rename({'time': 'time_counter','lat': 'y', 'lon': 'x'})
    # df2 = xr.open_dataset(ls2[idx])
    # df2 = df2.isel(lon=slice(448,938),lat=slice(176,321))
    # df2 = df2.reindex(lat=list(reversed(df2.lat))) 
    # df2 = df2.rename({'time': 'time_counter','lat': 'y', 'lon': 'x'})
    # # convert variables
    # # precip 
    # df['precip'] = df.precip*rho0/(pprec_fqh*3600)
    # df['precip'].attrs['units'] = 'kg/m2/s' 
    # df['SSR'] = df2['SSR']
