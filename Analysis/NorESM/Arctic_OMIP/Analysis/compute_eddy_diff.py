from mpl_toolkits.axes_grid1 import make_axes_locatable

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)



iind = 95
gr = xr.open_dataset('/cluster/shared/noresm/inputdata/ocn/micom/tnx1v4/20170601/grid.nc')

root_dir = '/tos-project1/NS2345K/noresm/cases/NOIIA_T62_tn11_sr10m60d_01/ocn/hist/'

yearstr = 281
yearend = 301
var = np.zeros((70,385))

for year in range(yearstr,yearend):
    for month in range(1,13):
        aa = f'{month:02}'
        fname = root_dir+'NOIIA_T62_tn11_sr10m60d_01.micom.hm.0'+np.str(year)+'-'+aa+'.nc'
        print(fname)
        df = xr.open_dataset(fname)['difisolvl']
        var = var + 10**df[0,:,:,iind]


var = var/240.0


root_dir = '/tos-project1/NS2345K/noresm/cases/NOIIAOC20TR_T62_tn14_20190628/ocn/hist/'

list = sorted(glob.glob(root_dir+'*hm*1988*'))
for year in range(1989,2008):
    list.extend(sorted(glob.glob(root_dir+'*hm*'+np.str(year)+'*')))


df = xr.open_mfdataset(list)['difisolvl']
var1 = df[:,:,:,iind]
var1 = var1.mean('time')


fig, (ax1, ax2) = plt.subplots(2)
# fig.suptitle('Vertically stacked subplots')
# ax1.plot(x, y)
# ax2.plot(x, -y)

img1 = ax1.pcolormesh(gr.plat[:200,iind],-var.depth,(var[:,:200].data),vmin=100,vmax=1500,cmap='needJet2',rasterized=True);
colorbar(img1)
# fig.colorbar(img1, ax=ax1)
# ax1.set_aspect('auto')
img2 = ax2.pcolormesh(gr.plat[:200,iind],-var.depth,(var1[:,:200].data),vmin=100,vmax=1500,cmap='needJet2',rasterized=True)
colorbar(img2)
# fig.colorbar(img2, ax=ax2)
# ax2.set_aspect('auto')

ax1.set(ylabel='Depth [m]');
ax2.set(ylabel='Depth [m]');
ax2.set(xlabel='Latitude');

# plt.savefig('paperfigs/Eddy_diff.eps', bbox_inches='tight',format='eps', dpi=200)
plt.savefig('paperfigs/Eddy_diff.png', bbox_inches='tight',format='png', dpi=300)
# plt.close()


