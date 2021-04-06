#!/usr/bin/env python
# coding: utf-8
import numpy as np
import xarray as xr
import pandas as pd
import xclim
import os
from netCDF4 import Dataset

# ------------------------------------------------------------------------------------------------------------------

# Data Path
input_path  = '/okyanus/users/etoker/CMIP6/Med_Data/day/'
output_path = '/okyanus/users/milicak/dataset/CMIP6/Python_Script/Clim_Indices/Output_Indice_nc/'

# Model Name

#1.  ACCESS-CM2
#2.  ACCESS-ESM1-5
#3.  BCC-CSM2-MR
#4.  CanESM5
#5.  EC-Earth3
#6.  EC-Earth3-Veg
#7.  FGOALS-g3      (variable: "tas" is missing on 2016/01/01 - 2019/12/31 and 2100/01/01 - 2100/12/31)
#8.  GFDL-ESM4
#9.  INM-CM4-8
#10. INM-CM5-0
#11. IPSL-CM6A-LR
#12. MIROC6
#13. MPI-ESM1-2-LR
#14. MRI-ESM2-0
#15. NESM3

# Ex: GFDL-ESM4 Mediterranean(Med) Daily Data
# model = 'GFDL-ESM4'

# all_models_name = ['GFDL-ESM4'] # you can write more than 1 model. Ex:  all_models_name = ['GFDL-ESM4','NESM3']
# all_models_name = ['CanESM5', 'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-g3',
all_models_name = ['GFDL-ESM4', 'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC6',
                   'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'ACCESS-CM2',
                   'ACCESS-ESM1-5', 'BCC-CSM2-MR']

for model in all_models_name:
    print(model)

    # Create specific output folder
    os.mkdir(output_path+model)

    # ------------------------------------------------------------------------------------------------------------------


    # Time Step

    # 251 Years x 365.4 days == 91676
    # 251 Years x 365 days == 91615 , without leap days (wld)

    #1.  ACCESS-CM2     - 91676
    #2.  ACCESS-ESM1-5  - 91676
    #3.  BCC-CSM2-MR    - 91615 - wld
    #4.  CanESM5        - 91615 - wld
    #5.  EC-Earth3      - 91676
    #6.  EC-Earth3-Veg  - 91676
    #7.  FGOALS-g3      - 91615 - wld , (variable: "tas" is missing on 2016/01/01 - 2019/12/31 and 2100/01/01 - 2100/12/31)
    #8.  GFDL-ESM4      - 91615 - wld
    #9.  INM-CM4-8      - 91615 - wld
    #10. INM-CM5-0      - 91615 - wld
    #11. IPSL-CM6A-LR   - 91676
    #12. MIROC6         - 91676
    #13. MPI-ESM1-2-LR  - 91676
    #14. MRI-ESM2-0     - 91676
    #15. NESM3          - 91676


    if model=='BCC-CSM2-MR' or model=='CanESM5' or model=='FGOALS-g3' or model=='GFDL-ESM4' or model=='INM-CM4-8' or model=='INM-CM5-0':
        time_step = pd.date_range('1850-01-01', freq='D', periods = 91676)
        time_step = time_step[ (time_step.month != 2) | (time_step.day != 29) ]    # delete leap day (29th day of 2nd month)
    else:
        time_step = pd.date_range('1850-01-01', freq='D', periods = 91676)



    # ------------------------------------------------------------------------------------------------------------------

    # Grid Str

    #1.  ACCESS-CM2     - gn
    #2.  ACCESS-ESM1-5  - gn
    #3.  BCC-CSM2-MR    - gn
    #4.  CanESM5        - gn
    #5.  EC-Earth3      - gr
    #6.  EC-Earth3-Veg  - gr
    #7.  FGOALS-g3      - gn  , (variable: "tas" is missing on 2016/01/01 - 2019/12/31 and 2100/01/01 - 2100/12/31)
    #8.  GFDL-ESM4      - gr1
    #9.  INM-CM4-8      - gr1
    #10. INM-CM5-0      - gr1
    #11. IPSL-CM6A-LR   - gr
    #12. MIROC6         - gn
    #13. MPI-ESM1-2-LR  - gn
    #14. MRI-ESM2-0     - gn
    #15. NESM3          - gn

    if model=='GFDL-ESM4' or model=='INM-CM4-8' or model=='INM-CM5-0':
        grid = 'gr1'
    elif model=='EC-Earth3' or model=='EC-Earth3-Veg' or model=='IPSL-CM6A-LR':
        grid = 'gr'
    else:
        grid = 'gn'

    # ------------------------------------------------------------------------------------------------------------------

    # Model Data

    pr_model_1850_2100_Med     = xr.open_dataset(input_path+'pr/'+model+'/pr_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')
    tas_model_1850_2100_Med    = xr.open_dataset(input_path+'tas/'+model+'/tas_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')
    tasmax_model_1850_2100_Med = xr.open_dataset(input_path+'tasmax/'+model+'/tasmax_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')
    tasmin_model_1850_2100_Med = xr.open_dataset(input_path+'tasmin/'+model+'/tasmin_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    pr_model_1850_2100_Med['time']     = time_step
    tas_model_1850_2100_Med['time']    = time_step
    tasmax_model_1850_2100_Med['time'] = time_step
    tasmin_model_1850_2100_Med['time'] = time_step

    p_model_1850_2100_Med  = pr_model_1850_2100_Med.pr
    t_model_1850_2100_Med  = tas_model_1850_2100_Med.tas
    tx_model_1850_2100_Med = tasmax_model_1850_2100_Med.tasmax
    tn_model_1850_2100_Med = tasmin_model_1850_2100_Med.tasmin

    # the base period
    p_model_1981_2010_Med  = p_model_1850_2100_Med.sel(time=slice('1981-01-01','2010-12-31'))
    t_model_1981_2010_Med  = t_model_1850_2100_Med.sel(time=slice('1981-01-01','2010-12-31'))
    tx_model_1981_2010_Med = tx_model_1850_2100_Med.sel(time=slice('1981-01-01','2010-12-31'))
    tn_model_1981_2010_Med = tn_model_1850_2100_Med.sel(time=slice('1981-01-01','2010-12-31'))


    # ------------------------------------------------------------------------------------------------------------------


    # ET-SCI indices - the sector(s)
    # H    =Health
    # AFS  =Agriculture and Food Security
    # WRH  =Water Resources and Hydrology

    # ETCCDI indices - http://etccdi.pacificclimate.org/list_27_indices.shtml


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 1 - FD - Frost Days (H, AFS)
    # Number of days where daily minimum temperatures are below 0℃. ( TN < 0 °C )
    # tasmin

    indice   = 'fd'
    longname = 'Frost Days'
    unit     = 'day'
    variable = 'tasmin'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.FD(tn_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.FD(tn_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 2 - SU - Summer Days (H)
    # Number of days where daily maximum temperature exceed a threshold. ( TX > 25 °C )
    # tasmax

    indice   = 'su'
    longname = 'Summer Days'
    unit     = 'day'
    variable = 'tasmax'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.SU(tx_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.SU(tx_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 3 - ID - Ice Days (H, AFS)
    # Number of ice/freezing days, where daily maximum temperatures are below 0℃. ( TX < 0 °C )
    # tasmax

    indice   = 'id'
    longname = 'Ice Days'
    unit     = 'day'
    variable = 'tasmax'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.ID(tx_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.ID(tx_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 4 - TR - Tropical Nights (H, AFS)
    # The number of days with minimum daily temperature above threshold. ( TN > 20 °C )
    # tasmin

    indice   = 'tr'
    longname = 'Tropical Nights'
    unit     = 'day'
    variable = 'tasmin'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TR(tn_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TR(tn_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 5 - GSL - Growing Season Length (AFS)
    # The number of days between the first occurrence of at least 6 consecutive days with mean daily
    # temperature over a threshold ( TM > 5 °C ) and the first occurrence of at least six consecutive days with mean
    # daily temperature below the same threshold ( TM < 5 °C ) after a certain date.
    # (Usually July 1st in the northern hemisphere and January 1st in the southern hemisphere.)
    # tas

    indice   = 'gsl'
    longname = 'Growing Season Length'
    unit     = 'day'
    variable = 'tas'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.GSL(t_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.GSL(t_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 6 - TXx - Max TX, Highest Maximum Temperature (AFS)
    # The maximum value of daily maximum temperature.
    # tasmax

    indice   = 'txx'
    longname = 'Highest Maximum Temperature'
    unit     = 'celsius'
    variable = 'tasmax'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TXx(tx_model_1850_2100_Med, freq = 'YS') - 273.15
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TXx(tx_model_1850_2100_Med, freq = 'MS') - 273.15
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 7 - TNn - Min TN, Lowest Minimum Temperature (AFS)
    # Minimum of daily minimum temperature.
    # tasmin

    indice   = 'tnn'
    longname = 'Lowest Minimum Temperature'
    unit     = 'celsius'
    variable = 'tasmin'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TNn(tn_model_1850_2100_Med, freq = 'YS') - 273.15
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TNn(tn_model_1850_2100_Med, freq = 'MS') - 273.15
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 8 - TXn - Min TX, Lowest Maximum Temperature
    # The minimum of daily maximum temperature.
    # tasmax

    indice   = 'txn'
    longname = 'Lowest Maximum Temperature'
    unit     = 'celsius'
    variable = 'tasmax'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TXn(tx_model_1850_2100_Med, freq = 'YS') - 273.15
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TXn(tx_model_1850_2100_Med, freq = 'MS') - 273.15
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 9 - TNx - Max TN, Highest Minimum Temperature
    # The maximum of daily minimum temperature.
    # tasmin

    indice   = 'tnx'
    longname = 'Highest Minimum Temperature'
    unit     = 'celsius'
    variable = 'tasmin'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TNx(tn_model_1850_2100_Med, freq = 'YS') - 273.15
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TNx(tn_model_1850_2100_Med, freq = 'MS') - 273.15
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 10 - TN10p - Amount of Cold Nights
    # Percentage of days with daily minimum temperature below the 10th percentile.
    # Let TN10 be the calendar day 10th percentile centred on a 5-day window for the base period 1981-2010
    # The percentage of time for the base period is determined where: TN < TN10
    # (day of year = 366, day of month = 31)
    # tasmin

    indice   = 'tn10p'
    longname = 'Amount of Cold Nights'
    unit     = '%'
    variable = 'tasmin'

    p = xclim.core.calendar.percentile_doy(tn_model_1981_2010_Med, window = 5, per=0.1)

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TN10p(tn_model_1850_2100_Med, p, freq = 'YS')
    indice_YS_model_1850_2100_Med = indice_YS_model_1850_2100_Med*100/365
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TN10p(tn_model_1850_2100_Med, p, freq = 'MS')
    indice_MS_model_1850_2100_Med = indice_MS_model_1850_2100_Med*100/31
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 11 - TX10p - Amount of Cool Days, Percentage of days when TX < 10th percentile
    # Number of days with daily maximum temperature below the 10th percentile.
    # Let TX10 be the calendar day 10th percentile centred on a 5-day window for the base period 1981-2010
    # The percentage of time for the base period is determined where: TX < TX10
    # tasmax

    indice   = 'tx10p'
    longname = 'Amount of Cool Days'
    unit     = '%'
    variable = 'tasmax'

    p = xclim.core.calendar.percentile_doy(tx_model_1981_2010_Med, window = 5, per=0.1)

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TX10p(tx_model_1850_2100_Med, p, freq = 'YS')
    indice_YS_model_1850_2100_Med = indice_YS_model_1850_2100_Med*100/365
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TX10p(tx_model_1850_2100_Med, p, freq = 'MS')
    indice_MS_model_1850_2100_Med = indice_MS_model_1850_2100_Med*100/31
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 12 - TN90p - Amount of Warm Nights, Percentage of days when TN > 90th percentile
    # Number of days with daily minimum temperature over the 90th percentile.
    # Let TN90 be the calendar day 90th percentile centred on a 5-day window for the base period 1981-2010
    # The percentage of time for the base period is determined where: TN > TN90
    # tasmin

    indice   = 'tn90p'
    longname = 'Amount of Warm Nights'
    unit     = '%'
    variable = 'tasmin'

    p = xclim.core.calendar.percentile_doy(tn_model_1981_2010_Med, window = 5, per=0.9)

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TN90p(tn_model_1850_2100_Med, p, freq = 'YS')
    indice_YS_model_1850_2100_Med = indice_YS_model_1850_2100_Med*100/365
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TN90p(tn_model_1850_2100_Med, p, freq = 'MS')
    indice_MS_model_1850_2100_Med = indice_MS_model_1850_2100_Med*100/31
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 13 - TX90p - Amount of Hot Days, Percentage of days when TX > 90th percentile
    # Number of days with daily maximum temperature over the 90th percentile.
    # Let TX90 be the calendar day 90th percentile centred on a 5-day window for the base period 1981-2010
    # The percentage of time for the base period is determined where: TX > TX90
    # tasmax

    indice   = 'tx90p'
    longname = 'Amount of Hot Days'
    unit     = '%'
    variable = 'tasmax'

    p = xclim.core.calendar.percentile_doy(tx_model_1981_2010_Med, window = 5, per=0.9)

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.TX90p(tx_model_1850_2100_Med, p, freq = 'YS')
    indice_YS_model_1850_2100_Med = indice_YS_model_1850_2100_Med*100/365
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.TX90p(tx_model_1850_2100_Med, p, freq = 'MS')
    indice_MS_model_1850_2100_Med = indice_MS_model_1850_2100_Med*100/31
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 14 - WSDI - Warm Speel Duration Index (H, AFS, WRH)
    # Number of days with at least six consecutive days where the daily maximum temperature is above the 90th percentile.
    # Let TX90 be the calendar day 90th percentile centred on a 5-day window for the base period 1981-2010
    # Then the number of days per period is summed where, in intervals of at least 6 consecutive days: TX > TX90
    # tasmax

    indice   = 'wsdi'
    longname = 'Warm Speel Duration Index'
    unit     = 'day'
    variable = 'tasmax'

    p = xclim.core.calendar.percentile_doy(tx_model_1981_2010_Med, window = 5, per=0.9)

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.WSDI(tx_model_1850_2100_Med, p, window = 6, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.WSDI(tx_model_1850_2100_Med, p, window = 6, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 15 - CSDI - Cold Speel Duration Index (H, AFS)
    # Number of days with at least six consecutive days where the daily minimum temperature is below the 10th percentile.
    # Let TN10 be the calendar day 10th percentile centred on a 5-day window for the base period 1981-2010.
    # Then the number of days per period is summed where, in intervals of at least 6 consecutive days: TN < TN10
    # tasmin

    indice   = 'csdi'
    longname = 'Cold Speel Duration Index'
    unit     = 'day'
    variable = 'tasmin'

    p = xclim.core.calendar.percentile_doy(tn_model_1981_2010_Med, window = 5, per=0.1)

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.CSDI(tn_model_1850_2100_Med, p, window = 6, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.CSDI(tn_model_1850_2100_Med, p, window = 6, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 16 - DTR - Daily Temperature Range
    # The mean difference between the daily maximum temperature and the daily minimum temperature.
    # tasmax,tasmin

    indice   = 'dtr'
    longname = 'Daily Temperature Range'
    unit     = 'celcius'
    variable = 'tasmax_tasmin'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.DTR(tx_model_1850_2100_Med, tn_model_1850_2100_Med, freq = 'YS') - 273.15
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.DTR(tx_model_1850_2100_Med, tn_model_1850_2100_Med, freq = 'MS') - 273.15
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 17 - RX1day - Monthly Maximum 1-day Precipitation
    # Resample the original daily total precipitation series by taking the max over each period.
    # pr

    indice   = 'rx1day'
    longname = 'Monthly Maximum 1-day Precipitation'
    unit     = 'mm'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.RX1day(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.RX1day(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 18 - RX5day - Monthly Maximum 5-day Precipitation
    # Calculate the 5-day rolling sum of the original daily total precipitation series and determine the maximum value over each period.
    # pr

    indice   = 'rx5day'
    longname = 'Monthly Maximum 5-day Precipitation'
    unit     = 'mm'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.RX5day(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.RX5day(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 19 - SDII - Simple Pricipitation Intensity Index
    # Return the average precipitation over wet days (PR ≥ 1mm)
    # pr

    indice   = 'sdii'
    longname = 'Simple Pricipitation Intensity Index'
    unit     = 'mm/day'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.SDII(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.SDII(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 20 - R10mm - Number of Heavy Rain Days
    # Return the total number of days during period with precipitation over threshold (PR ≥ 10mm)
    # pr

    indice   = 'r10mm'
    longname = 'Number of Heavy Rain Days'
    unit     = 'day'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.R10mm(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.R10mm(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 21 - R20mm - Number of Very Heavy Rain Days (AFS, WRH)
    # Return the total number of days during period with precipitation over threshold ( PR ≥ 20 mm )
    # pr

    indice   = 'r20mm'
    longname = 'Number of Heavy Rain Days'
    unit     = 'day'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.R20mm(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.R20mm(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 23 - CDD - Consecutive Dry Days (H, AFS, WRH)
    # Maximum number of consecutive dry days, longest dry spell
    # Return the maximum number of consecutive days within the period where
    # precipitation is below a certain threshold ( PR < 1.0 mm )
    # pr

    indice   = 'cdd'
    longname = 'Consecutive Dry Days'
    unit     = 'day'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.CDD(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.CDD(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 24 - CWD - Consecutive Wet Days (H, AFS, WRH)
    # Maximum number of consecutive wet days, longest wet spell
    # Return the maximum number of consecutive days within the period where
    # precipitation is over a certain threshold ( PR ≥ 1.0 mm )
    # pr

    indice   = 'cwd'
    longname = 'Consecutive Wet Days'
    unit     = 'day'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.CWD(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.CWD(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 25 - R95pTOT - Contribution from Very Wet Days (AFS, WRH)
    # Percentage of the total precipitation over period occurring in days
    # where the precipitation is above a threshold ( PR ≥ 1.0 mm ) defining wet days
    # and above a given percentile ( PR ≥ PR95 ) for that day
    # pr

    indice   = 'r95ptot'
    longname = 'Contribution from Very Wet Days'
    unit     = '%'
    variable = 'pr'

    p_model_1850_2100_Med_bt1mm = p_model_1850_2100_Med.where((p_model_1850_2100_Med*86400)>=1) # ( PR ≥ 1.0 mm )
    p = xclim.core.calendar.percentile_doy(p_model_1850_2100_Med_bt1mm, window = 1, per=0.95) # ( PR ≥ PR95 )

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.R95pTOT(p_model_1850_2100_Med, per=p, freq = 'YS')*100
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.R95pTOT(p_model_1850_2100_Med, per=p, freq = 'MS')*100
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------


    # Indice 26 - R99pTOT - Contribution from Extremely Wet Days (AFS, WRH)
    # Percentage of the total precipitation over period occurring in days
    # where the precipitation is above a threshold ( PR ≥ 1.0 mm ) defining wet days
    # and above a given percentile ( PR ≥ PR99 ) for that day
    # pr

    indice   = 'r99ptot'
    longname = 'Contribution from Extremely Wet Days'
    unit     = '%'
    variable = 'pr'

    p_model_1850_2100_Med_bt1mm = p_model_1850_2100_Med.where((p_model_1850_2100_Med*86400)>=1) # ( PR ≥ 1.0 mm )
    p = xclim.core.calendar.percentile_doy(p_model_1850_2100_Med_bt1mm, window = 1, per=0.99) # ( PR ≥ PR99 )

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.R99pTOT(p_model_1850_2100_Med, per=p, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.R99pTOT(p_model_1850_2100_Med, per=p, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')



    # ------------------------------------------------------------------------------------------------------------------


    # Indice 27 - PRCPTOT - Annual Total Wet-Day (AFS, WRH)
    # Accumulated total (liquid and(default)/or solid) precipitation. ( PR ≥ 1.0 mm )
    # (When the mean temperature is over 0 degC, precipitatio is assumed to be liquid rain and snow otherwise. Optionally, you can add temp data.)
    # pr

    indice   = 'prcptot'
    longname = 'Annual Total Wet-Day'
    unit     = 'mm'
    variable = 'pr'

    # Yearly
    indice_YS_model_1850_2100_Med = xclim.ICCLIM.PRCPTOT(p_model_1850_2100_Med, freq = 'YS')
    indice_YS_model_1850_2100_Med.name                   = indice
    indice_YS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_YS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_YS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_YS_model_1850_2100_Med.attrs['units']         = unit
    indice_YS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_YS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')

    # Monthly
    indice_MS_model_1850_2100_Med = xclim.ICCLIM.PRCPTOT(p_model_1850_2100_Med, freq = 'MS')
    indice_MS_model_1850_2100_Med.name                   = indice
    indice_MS_model_1850_2100_Med.attrs['original_name'] = indice
    indice_MS_model_1850_2100_Med.attrs['long_name']     = longname
    indice_MS_model_1850_2100_Med.attrs['standard_name'] = longname
    indice_MS_model_1850_2100_Med.attrs['units']         = unit
    indice_MS_model_1850_2100_Med.to_netcdf(output_path+model+'/'+indice+'_xclim_MS_'+variable+'_day_'+model+'_hist_ssp585_'+grid+'_185001-210012_Med.nc')


    # ------------------------------------------------------------------------------------------------------------------




