import xarray as xr
import matplotlib.pyplot as plt

topo=xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/2X/ocean_topog.nc')
Omask=xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/2X/ocean_mask.nc')
Lmask=xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/2X/land_mask.nc')
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

# select and ocean point
new_mask=ice9(100,100,Omask.mask.data,False,False)
newOmask=xr.Dataset()
newLmask=xr.Dataset()
newOmask['mask']=xr.DataArray(new_mask, coords=Omask.coords)
newLmask['mask']=xr.DataArray(new_mask, coords=Omask.coords)
print(newOmask)

newLmask['mask']=xr.where(newOmask.mask==0.0,1.0,0.0)


inputDir='./'

topog2=topo.assign_coords(ntiles=('ntiles',[1]))
aa=newOmask.mask
aa=aa.rename({'nx':'ny','ny':'nx'})
topog2['depth']=topog2.depth.where(aa>0)
topog2['depth']=topog2.depth.fillna(0)
topog2.to_netcdf(f'{inputDir}ocean_topog.nc', mode='w', format='NETCDF3_64BIT')
newLmask.to_netcdf(f'{inputDir}land_mask.nc', mode='w', format='NETCDF3_64BIT')
newOmask.to_netcdf(f'{inputDir}ocean_mask.nc', mode='w', format='NETCDF3_64BIT')


