import numpy as np
import xarray as xr
import scipy
from datetime import datetime, timedelta
from my_sea_over_land import seaoverland_2d,seaoverland_1d 
from shyfem_utils import write_record_in_file1d, write_record_in_file2d, write_record_in_file2d_vel


root_folder = '/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/analysis/day/'
date0 = datetime(2020,1,1,0)
# dates = [date0 + timedelta(days=d) for d in range(366)]
dates = [date0 + timedelta(days=d) for d in range(31)]
ssh_shift = 0.0 #0.3


date_string = '%Y%m%d'
# root_filename = '%s/%s/%s_d-CMCC--%s-BSeas3-BS-b20190101_an-sv09.00.nc'
# root_filename = '%s/%s/%s_d-CMCC--%s-BSeas3-BS-b2019*_an-fv07.00.nc'
root_filename = '%s/%s/%s_d-CMCC--%s-BSeas4-BS-b2021*_an-sv10.00.nc'
outnetcdfname = 'uTSS_boundary_valuesL3.nc'

boundary_nodes_path = 'uTSS_boundary_list_lonlat3.txt'
out_file_temp = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/tempn_L3.dat'
varname_temp  = ' temperature [C]\n'
out_file_salt = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/saltn_L3.dat'
varname_salt  = ' salinity [psu]\n'
out_file_ssh  = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/boundn_L3.dat'
varname_ssh   = ' water level [m]\n'
out_file_uv   = '/work/opa/mi19918/Projects/uTSS_SHYFEM/forcing_files/uv3d_L3.dat'
varname_uvel  = ' u-velocity [m/s]\n'
varname_vvel  = ' v-velocity [m/s]\n'

data = np.genfromtxt(boundary_nodes_path,usecols=(0,1))
boundary_lons,boundary_lats = data[:,0],data[:,1]

n_bound = len(boundary_lons)

## create target interpolation depths
utss_bottom_levels= [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 
                    16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 
                    28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 
                    41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 55.0, 60.0, 65.0, 
                    70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 120.0, 140.0, 160.0,
               180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0,
               1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0 ]



utss_level_interfaces = np.roll(np.append(np.asarray(utss_bottom_levels),0),1)
interpolation_levels = 0.5*(utss_level_interfaces[1:] + utss_level_interfaces[:-1])

n_lev = len(utss_bottom_levels)
x = xr.DataArray(boundary_lons,dims='boundary_nodes')
y = xr.DataArray(boundary_lats,dims='boundary_nodes')


fout_temp = open(out_file_temp,'w')
fout_salt = open(out_file_salt,'w')
fout_ssh  = open(out_file_ssh,'w')
fout_uv   = open(out_file_uv,'w')

container = xr.Dataset()

for date in dates:
    print ('loading date %s' % date.strftime('%Y/%m/%d'))
    if date.strftime('%m')=='12' and date.strftime('%d')=='23':
        root_filename = '%s/%s/%s_d-CMCC--%s-BSeas3-BS-b2020*_an-sv09.00.nc'
    
    root_filename2 = root_folder + root_filename
    temp_filename = root_filename2 % (date.strftime('%Y'),date.strftime('%m'),date.strftime(date_string),'TEMP')
    salt_filename = root_filename2 % (date.strftime('%Y'),date.strftime('%m'),date.strftime(date_string),'PSAL')
    vel_filename  = root_filename2 % (date.strftime('%Y'),date.strftime('%m'),date.strftime(date_string),'RFVL')
    ssh_filename  = root_filename2 % (date.strftime('%Y'),date.strftime('%m'),date.strftime(date_string),'ASLV')
    
    dnm = sorted(glob.glob(temp_filename))
    temp_filename = dnm[0]
    dnm = sorted(glob.glob(salt_filename))
    salt_filename = dnm[0]
    dnm = sorted(glob.glob(vel_filename))
    vel_filename = dnm[0]
    dnm = sorted(glob.glob(ssh_filename))
    ssh_filename = dnm[0]
    
    ds_temp  = xr.open_dataset(temp_filename)
    ds_salt  = xr.open_dataset(salt_filename)
    ds_vel   = xr.open_dataset(vel_filename)
    ds_ssh   = xr.open_dataset(ssh_filename)
    
    ds_temp2 = ds_temp.interp(lon=x,lat=y,depth=interpolation_levels)
    ds_salt2 = ds_salt.interp(lon=x,lat=y,depth=interpolation_levels)
    ds_vel2  = ds_vel.interp(lon=x,lat=y,depth=interpolation_levels)
    ds_ssh2  = ds_ssh.interp(lon=x,lat=y)
    
    # do sea over land on transect
    mask2d = np.isnan(ds_salt2.so)
    mask1d = np.isnan(ds_ssh2.zos)
    transect_temp = np.ma.masked_array(ds_temp2.thetao.data,mask=mask2d)
    transect_salt = np.ma.masked_array(ds_salt2.so.data,mask=mask2d)
    transect_uvel = np.ma.masked_array(ds_vel2.uo.data,mask=mask2d)
    transect_vvel = np.ma.masked_array(ds_vel2.vo.data,mask=mask2d)
    transect_ssh  = np.ma.masked_array(ds_ssh2.zos.data,mask=mask1d)
    
    transect_salt_sol  = seaoverland_2d(transect_salt[0,:],20,copy=True)
    transect_temp_sol  = seaoverland_2d(transect_temp[0,:],20,copy=True)
    transect_uvel_sol  = seaoverland_2d(transect_uvel[0,:],20,copy=True)
    transect_vvel_sol  = seaoverland_2d(transect_vvel[0,:],20,copy=True)
    transect_ssh_sol   = seaoverland_1d(transect_ssh[0,:],10,copy=True)
    
    ## assign new variables with Sea Over Land to datasets, besides original transects
    #ds_temp2 = ds_temp2.assign(thetao_sol=ds_temp2["thetao"])
    #ds_salt2 = ds_salt2.assign(so_sol=ds_salt2["so"])
    #ds_vel2  = ds_vel2.assign(uo_sol=ds_vel2["uo"])
    #ds_vel2  = ds_vel2.assign(vo_sol=ds_vel2["vo"])
    #ds_ssh2  = ds_ssh2.assign(zos_sol=ds_ssh2["zos"])
    
    #ds_temp2.thetao_sol.data[0,:] = transect_temp_sol
    #ds_salt2.so_sol.data[0,:]     = transect_salt_sol
    #ds_vel2.uo_sol.data[0,:]      = transect_uvel_sol
    #ds_vel2.vo_sol.data[0,:]      = transect_vvel_sol
    #ds_ssh2.zos_sol.data[0,:]     = transect_ssh_sol
    
    ds_temp2.thetao.data[0,:] = transect_temp_sol
    ds_salt2.so.data[0,:]     = transect_salt_sol
    ds_vel2.uo.data[0,:]      = transect_uvel_sol
    ds_vel2.vo.data[0,:]      = transect_vvel_sol
    ds_ssh2.zos.data[0,:]     = transect_ssh_sol
    
    ## write recod in files
    write_record_in_file1d(fout_ssh,date,varname_ssh,n_bound,1,1,np.copy(transect_ssh_sol)+ssh_shift)
    write_record_in_file2d(fout_temp,date,utss_bottom_levels,varname_temp,n_bound,n_lev,1,np.copy(transect_temp_sol))
    write_record_in_file2d(fout_salt,date,utss_bottom_levels,varname_salt,n_bound,n_lev,1,np.copy(transect_salt_sol))
    write_record_in_file2d_vel(fout_uv,date,utss_bottom_levels,varname_uvel,varname_vvel,n_bound,n_lev,2,
                               np.copy(transect_uvel_sol),np.copy(transect_vvel_sol))
    
    ds_merged = xr.merge([ ds_temp2, ds_salt2, ds_vel2, ds_ssh2  ])
    
    if date == date0:
    	container = ds_merged
    else:
    	container = xr.concat([container,ds_merged],dim='time')
    


container.to_netcdf(outnetcdfname)

fout_temp.close()
fout_salt.close()
fout_ssh.close()
fout_uv.close()

