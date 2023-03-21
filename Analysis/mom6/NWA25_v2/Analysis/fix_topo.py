import xarray as xr
import matplotlib.pyplot as plt

topo=xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_new/ocean_topog.nc')
Omask=xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_new/ocean_mask.nc')
Lmask=xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_new/land_mask.nc')
topo.depth.plot(vmin=0)

#https://github.com/raphaeldussin/om4labs/blob/b62fb4ca85516997911f024b236087b61d65a039/om4labs/m6toolbox.py#L173
def ice9(i, j, source, xcyclic=True, tripolar=True):
    """
    An iterative (stack based) implementation of "Ice 9".
    The flood fill starts at [j,i] and treats any positive value of "source" as
    passable. Zero and negative values block flooding.
    xcyclic = True allows cyclic behavior in the last index. (default)
    tripolar = True allows a fold across the top-most edge. (default)
    Returns an array of 0's and 1's.
    """
    wetMask = 0 * source
    (nj, ni) = wetMask.shape
    stack = set()
    stack.add((j, i))
    while stack:
        (j, i) = stack.pop()
        if wetMask[j, i] or source[j, i] <= 0:
            continue
        wetMask[j, i] = 1
        if i > 0:
            stack.add((j, i - 1))
        elif xcyclic:
            stack.add((j, ni - 1))
        if i < ni - 1:
            stack.add((j, i + 1))
        elif xcyclic:
            stack.add((j, 0))
        if j > 0:
            stack.add((j - 1, i))
        if j < nj - 1:
            stack.add((j + 1, i))
        elif tripolar:
            stack.add((j, ni - 1 - i))  # Tri-polar fold
    return wetMask


new_mask=ice9(500,200,Omask.mask.data,False,False)
newOmask=xr.Dataset()
newLmask=xr.Dataset()
newOmask['mask']=xr.DataArray(new_mask, coords=Omask.coords)
newLmask['mask']=xr.DataArray(new_mask, coords=Omask.coords)
print(newOmask)


#Corrections

newLmask['mask']=xr.where(newOmask.mask==0.0,1.0,0.0)


inputDir='./'

topog2=topo.assign_coords(ntiles=('ntiles',[1]))
topog2['depth']=topog2.depth.where(newOmask.mask>0)
topog2['depth']=topog2.depth.fillna(0)
topog2.to_netcdf(f'{inputDir}ocean_topog.nc', mode='w', format='NETCDF3_64BIT')
newLmask.to_netcdf(f'{inputDir}land_mask.nc', mode='w', format='NETCDF3_64BIT')
newOmask.to_netcdf(f'{inputDir}ocean_mask.nc', mode='w', format='NETCDF3_64BIT')


from gridtools.gridutils import GridUtils
import sys, os, logging, cartopy

inputDir = "/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_new/"

# Initialize a grid object
grd = GridUtils()
grd.openGrid(f"{inputDir}ocean_hgrid.nc", gridType='MOM6')
grd.readGrid()

# Write out FMS related support files
grd.makeSoloMosaic(
    topographyGrid=topo['depth'],
    writeLandmask=True,
    writeOceanmask=True,
    inputDirectory=inputDir,
    overwrite=True,
)
grd.saveGrid(filename=os.path.join(inputDir, "ocean_hgrid_final.nc"))


topo_check=xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_new/ocean_topog.nc')
p=topo_check.depth.plot( size=10, cmap= 'gist_earth', levels=np.arange(0,8000,1000))
plt.contour(topo_check.depth,colors='k',levels=np.arange(5000,8000,1000))


x = 56; y = 200
dx=50
topo_check.depth.isel(nx=slice(x-dx,x+dx+1),ny=slice(y-dx,y+dx+1)).plot(
    cmap= 'gist_earth', levels=np.arange(0,6000,200), size=7)
plt.contour(topo_check.depth.isel(nx=slice(x-dx,x+dx+1),ny=slice(y-dx,y+dx+1))
                                  ,colors='k',levels=[50,100,150,200,600,1200,2400,3600,4800,5000])


dx=15
topo_check.depth.isel(nx=slice(x-dx,x+dx+1),ny=slice(y-dx,y+dx+1)).plot(
    cmap= 'gist_earth', levels=np.arange(50,170,5), size=7)
plt.contour(topo_check.depth.isel(nx=slice(x-dx,x+dx+1),ny=slice(y-dx,y+dx+1))
                                  ,colors='k',levels=np.arange(50,150,5))


