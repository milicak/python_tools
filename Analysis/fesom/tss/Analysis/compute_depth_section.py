import numpy as np
import sys

sys.path.append('/okyanus/users/milicak/python_libs/pyfesom')
import pyfesom as pf

meshpath = '/okyanus/users/milicak/models/fesom2/FESOM2_minimum_input/mesh/tss_Ali'

mesh = pf.load_mesh(meshpath, abg=[0, 0, 0], usepickle=True)
lon = mesh.x2
lat = mesh.y2

df = xr.open_dataset('/okyanus/users/milicak/Analysis/fesom/TSS/Analysis/oper.2018.oce.mean.nc')

level_data, elem_no_nan = pf.get_data(df['v'][0,:], mesh, 750)

plt.tripcolor(lon,lat,elem_no_nan[::], level_data,cmap='needJet2');plt.colorbar();

