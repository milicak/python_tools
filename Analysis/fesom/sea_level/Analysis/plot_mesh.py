import sys
sys.path.append('/okyanus/users/milicak/python_libs/pyfesom')
import pyfesom as pf
from mpl_toolkits.basemap import Basemap


meshpath='/okyanus/users/milicak/models/fesom2/FESOM2_minimum_input/mesh/mesh_tss_final/'
mesh = pf.load_mesh(meshpath, get3d=False, abg=[0,0,0])
elem2=mesh.elem[mesh.no_cyclic_elem,:]
m = Basemap(projection='robin',lon_0=0, resolution='l')
x, y = m(mesh.x2, mesh.y2)


plt.figure(figsize=(9,5))
m.drawmapboundary(fill_color='0.9')
m.drawcoastlines()
plt.triplot(x, y, elem2, lw=0.2,color='k');

