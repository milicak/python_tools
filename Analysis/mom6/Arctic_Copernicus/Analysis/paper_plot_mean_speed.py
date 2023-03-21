import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
from xgcm import Grid
from pandas.tseries.offsets import DateOffset

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'
deg_rad = np.pi/180.
ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[0:12]
for itr in range(2,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2[24:],decode_times=False)
time = pd.date_range("1998-01-15", freq=DateOffset(months=1), periods=ds.Time.shape[0])
df['Time'] = time

gr = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_geometry.nc')
lon = np.copy(gr.geolon)
lat = np.copy(gr.geolat)
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})

gr2 = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_hgrid.nc')
angle_dx = gr2.angle_dx[1::2,1::2]

# merge data
df = df.merge(gr)

# 2D grid
grid = Grid(df, coords={'X': {'center': 'xh', 'outer': 'xq'},
                        'Y': {'center': 'yh', 'outer': 'yq'},
                         }, periodic=[])

urho = grid.interp(df.SSU, 'X', boundary='fill')
vrho = grid.interp(df.SSV, 'Y', boundary='fill')

umean = urho.mean('Time')
vmean = vrho.mean('Time')

ue = np.cos(angle_dx.data*deg_rad)*umean.data+np.sin(angle_dx.data*deg_rad)*vmean.data
vn = -np.sin(angle_dx.data*deg_rad)*umean.data+np.cos(angle_dx.data*deg_rad)*vmean.data

sp = np.sqrt(ue**2+vn**2)

m = Basemap(projection='npstere',boundinglat=42,lon_0=0,resolution='i');
lon1,lat1 = m(-132,26)

fig, axes = plt.subplots(figsize=(6,6))
ax1 = plt.subplot2grid(shape=(1,1),loc=(0,0),colspan=1)
plt.tight_layout()

m.ax=ax1
m.drawcoastlines()
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(gr.D==0,sp),
                   cmap='Blues_r',vmin=0,vmax=0.3,latlon=True, rasterized=True)
# parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
# m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
# meridians = np.arange(10.,351.,20.)
# m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.05,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im1, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r'm/s',rotation=0,y=1.06,labelpad=-45,fontsize=14)
# ax1.text(lon1,lat1,'a)',fontsize=14)

plt.savefig('paperfigs/mean_speed.png', bbox_inches='tight',format='png',dpi=300)

