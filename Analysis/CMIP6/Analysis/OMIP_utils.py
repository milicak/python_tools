import xarray as xr
import numpy as np
import os
import glob
#from joblib import Parallel, delayed
from dask import compute, delayed

def upper_700m_mean(templvl, pbot, pdepth, pmask, depth_bnds):
    ''' '''
    #
    dum700      = xr.ones_like(templvl.isel(time=0)).drop('time').rename('dum700')*700
    depth_bnds1 = xr.ones_like(dum700)*depth_bnds.isel(bounds=1).isel(time=0).drop('time')
    depth_bnds0 = xr.ones_like(dum700)*depth_bnds.isel(bounds=0).isel(time=0).drop('time')
    pdepth3D    = xr.ones_like(dum700)*pdepth
    #
    dplvl       = (xr.concat([depth_bnds1,dum700,pdepth3D],dim='dum').min(dim='dum') - xr.concat([depth_bnds0,dum700,pdepth3D],dim='dum').min(dim='dum'))*pbot/pdepth 
    #
    thetao700   = (templvl*dplvl*pmask).sum(dim='depth')/(dplvl*pmask).sum(dim='depth')
    #
    return thetao700.rename('thetao')

def thermosteric_sealvl_ref_period(th,s,dz,dp,parea,pmask,reference_slice=slice(0,5*12)):
    '''
    Calculate annual mean thermosteric sealevel component
    
    the assumption is that the reference period is within the full period 
    '''
    #
    V0   = annual_mean((dz.isel(time=reference_slice)*parea*pmask).sum(dim=('sigma','x','y'))).mean(dim='time').rename('V0')
    A    = (pmask*parea).sum().rename('A')
    #
    p0   = dp.isel(time=reference_slice).cumsum(dim='sigma')
    s0   = s.isel(time=reference_slice)
    th0  = th.isel(time=reference_slice)
    rho0 = eosben07_rho(p0,s0,th0)
    rho0 = annual_mean((parea*dp.isel(time=reference_slice)*rho0).sum(dim=('sigma','x','y'))/(parea*dp.isel(time=reference_slice)).sum(dim=('sigma','x','y'))).mean(dim='time').rename('rho0')
    p0   = annual_mean(p0).mean(dim='time').rename('p0')
    s0   = annual_mean(s0).mean(dim='time').rename('s0')
    th0  = annual_mean(th0).mean(dim='time').rename('th0')
    xr.merge([p0,s0,th0,rho0,V0,A]).to_netcdf('reference_period.nc')
    ref_per = xr.open_dataset('reference_period.nc')
    #
    #rho  = eosben07_rho(ref_per.p0,ref_per.s0,th)
    #rho  = annual_mean((parea*dp*rho).sum(dim=('sigma','x','y'))/(parea*dp).sum(dim=('sigma','x','y')))
    #
    #zostoga = (ref_per.V0/ref_per.A)*(1-rho/ref_per.rho0)
    #
    #return zostoga.rename('zostoga')
    return ref_per

def thermosteric_sealvl(th,dp,parea,ref_per):
    #
    rho  = eosben07_rho(ref_per.p0,ref_per.s0,th)
    rho  = (parea*dp*rho).sum(dim=('sigma','x','y'))/(parea*dp).sum(dim=('sigma','x','y'))
    #
    zostoga = (ref_per.V0/ref_per.A)*(1-rho/ref_per.rho0)
    #
    return zostoga.rename('zostoga')

def annual_mean(var,weights=np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])/365):
    '''
    Calculate annual means from monthly means assuming no-leap calendar
    '''
    month_weights = xr.DataArray(np.tile(weights,len(var.time)//12),coords=[var.time], name='month_weights')
    annual_mean = (month_weights*var).groupby('time.year').sum('time')
    #
    return annual_mean.rename({'year':'time'})

def annual_mean_to_file(var,fname,weights=np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])/365):
    '''
    Calculate annual means from monthly means assuming no-leap calendar
    '''
    month_weights = xr.DataArray(np.tile(weights,len(var.time)//12),coords=[var.time], name='month_weights')
    annual_mean = (month_weights*var).groupby('time.year').sum('time')
    annual_mean = annual_mean.rename({'year':'time'})
    annual_mean.rename(var.name).to_dataset().to_netcdf(fname)

def annual_mean_loop(var,dum_folder='/cluster/work/users/cgu025/tmp/',n_jobs=6,savepath=None):
    '''
    Calculates annual mean year-by-year, saves the individual years to a temprorary file, and opens all the files in one variable
    '''
    #
    if not os.path.exists(dum_folder):
        os.makedirs(dum_folder)
    elif len(glob.glob(dum_folder+'*.nc'))>0:
        for remfile in glob.glob(dum_folder+'*.nc'):
            os.remove(remfile)
    # calculate the annual mean year by year to save memory - then load the whole thing back - at this point no computations should be needed
    #Parallel(n_jobs=n_jobs)(delayed(annual_mean_to_file)(var.isel(time=slice(j*12,(j+1)*12)),dum_folder+var.name+'_dum_'+str(j).zfill(4)+'.nc') for j in range(len(var.time)//12))
    #values  = [delayed(annual_mean_to_file)(var.isel(time=slice(j*12,(j+1)*12)),dum_folder+var.name+'_dum_'+str(j).zfill(4)+'.nc') for j in range(len(var.time)//12)]
    #results = compute(*values)
    for j in range(len(var.time)//12):
        print(j)
        annual_mean_to_file(var.isel(time=slice(j*12,(j+1)*12)),dum_folder+var.name+'_dum_'+str(j).zfill(4)+'.nc')
    #
    var_out = xr.open_mfdataset(dum_folder+var.name+'_dum_*.nc',parallel=True)
    var_out=var_out.where(var_out!=0)  # cguo
    if savepath!=None:
        var_out.to_netcdf(savepath)
    #
    exec('var_out = var_out.'+var.name)
    #
    return var_out

def eosben07_const():
    a11= 9.9985372432159340e+02
    a12= 1.0380621928183473e+01
    a13= 1.7073577195684715e+00  
    a14=-3.6570490496333680e-02
    a15=-7.3677944503527477e-03
    a16=-3.5529175999643348e-03
    a21= 1.0  
    a22= 1.0316374535350838e-02
    a23= 8.9521792365142522e-04
    a24=-2.8438341552142710e-05
    a25=-1.1887778959461776e-05
    a26=-4.0163964812921489e-06
    b11= 1.7083494994335439e-02
    b12= 7.1567921402953455e-05
    b13= 1.2821026080049485e-05
    b21= 1.1995545126831476e-05
    b22= 5.5234008384648383e-08
    b23= 8.4310335919950873e-09
    
    return a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23

def p_alpha(p,p0,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    b1i=1/(b11+b12*th+b13*s)
    a1=(a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s))*b1i
    a2=(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s))*b1i
    b2=(b21+b22*th+b23*s)*b1i
    #
    r=b2*(p-p0)+(a2-a1*b2)*np.log((a1+p)/(a1+p0))
    #
    return r

def eosben07_rho(p,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    r=(a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s)+p*(b11+b12*th+b13*s))/(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s)+p*(b21+b22*th+b23*s))
    #
    return r

def eosben07_rho_th(p,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    P=a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s)+p*(b11+b12*th+b13*s)
    Qi=1/(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s)+p*(b21+b22*th+b23*s))
    #
    r=Qi*(a12+2.0*a14*th+a15*s+b12*p-Qi*P*(a22+2.0*a24*th+a25*s+b22*p))
    #
    return r

def eosben07_rho_s(p,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    P=a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s)+p*(b11+b12*th+b13*s);
    Qi=1/(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s)+p*(b21+b22*th+b23*s))
    #
    r=Qi*(a13+a15*th+2.0*a16*s+b13*p-Qi*P*(a23+a25*th+2.0*a26*s+b23*p))
    #
    return r
