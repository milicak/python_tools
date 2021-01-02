import gsw
from HCtFlood.kara import flood_kara
from scipy.signal import convolve as scipy_convolve
from scipy.signal import convolve2d as scipy_convolve2d


def strideConv(arr, arr2, s):
    return signal_convolve2d(arr, arr2[::-1, ::-1], mode='valid')[::s, ::s]

def mask_var(df,varname):
    ds = df[varname]
    mask = ds.where(ds!=0)
    mask = mask.where(np.isnan(mask)==1,1)
    df[varname] = ds*mask
    return df

# kernel 4 points
k = np.array([[0,1,0],[1,2,1],[0,1,0]])
# kernel to 1 degree
k = np.ones((10,10))
# kernel grad in x-direction Prewitt
kgrad_x = np.array([[-1,0,1],[-1,0,1],[-1,0,1]])
# kernel grad in y-direction Prewitt
kgrad_y = np.array([[1,1,1],[0,0,0],[-1,-1,-1]])


k = k/np.sum(k)

root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp01_0'
print(project_name)
outname = project_name + '_delta_rho.nc'
fnames = root_folder+project_name+'/*SSS*'
list = sorted(glob.glob(fnames))

# grname = root_folder+project_name+'/grid.nc'

# dfs = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/SSS_2010_08_14.nc')
# dft = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/SST_2010_08_14.nc')

# ds = flood_kara(dfs['SALT'], xdim='i', ydim='j', zdim='k')

# plt.pcolormesh(rhomean1-rhomean2,cmap='RdBu_r');plt.colorbar();

# mean new temperature, salinity
comp = dict(zlib=True, complevel=5)

for ind in np.arange(0,len(list)):
    print(ind)
    inname = list[ind][-13:]
    fname = root_folder + project_name + '/' + 'SST_' + inname
    dft = xr.open_dataset(fname)
    fname = root_folder + project_name + '/' + 'SSS_' + inname
    dfs = xr.open_dataset(fname)
    # mask salinity
    dfs =  mask_var(dfs,'SALT')
    # mask temperature
    dft =  mask_var(dft,'THETA')
    rho1 = gsw.sigma0(dfs.SALT,dft.THETA)
    # DO not forget to flip the kernel
    Tmean = scipy_convolve(dft.THETA[0,:,:],k[::-1, ::-1],mode='same', method='direct')
    Smean = scipy_convolve(dfs.SALT[0,:,:],k[::-1, ::-1],mode='same', method='direct')
    rhomean1 = scipy_convolve(rho1[0,:,:],k[::-1,::-1],mode='same', method='direct')
    rhomean2 = gsw.sigma0(Smean,Tmean)
    # delta_rho
    delta_rho = rhomean1-rhomean2
    Tmean_x = scipy_convolve(Tmean,kgrad_x[::-1,::-1],mode='same', method='direct')
    Tmean_y = scipy_convolve(Tmean,kgrad_y[::-1,::-1],mode='same', method='direct')
    Smean_x = scipy_convolve(Smean,kgrad_x[::-1,::-1],mode='same', method='direct')
    Smean_y = scipy_convolve(Smean,kgrad_y[::-1,::-1],mode='same', method='direct')
    # creata arrays
    drho = xr.DataArray(delta_rho, dims=("i", "j"))
    S_x = xr.DataArray(Smean_x, dims=("i", "j"))
    S_y = xr.DataArray(Smean_y, dims=("i", "j"))
    T_x = xr.DataArray(Tmean_x, dims=("i", "j"))
    T_y = xr.DataArray(Tmean_y, dims=("i", "j"))
    ds = xr.Dataset({"delta_rho": drho, "Salt_x": S_x , "Salt_y": S_y,
                     "Temp_x": T_x, "Temp_y": T_y})
    fname = root_folder + project_name + '/' + 'Delta_rho_' + inname
    encoding = {var: comp for var in ds.data_vars}
    print(fname)
    ds.to_netcdf(fname, encoding=encoding)



