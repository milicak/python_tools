import numpy as np


ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/2014*daily*'))

df = xr.open_mfdataset(ls1)
root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'

gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')


itr = 1
fig = plt.figure(figsize=(16,8))
m = Basemap(projection='merc',llcrnrlat=5,urcrnrlat=55,
            llcrnrlon=-100,urcrnrlon=-35,lat_ts=20,resolution='i')
for ind in range(90,360,4):
    print(itr)
    ax = fig.add_subplot(121)
    m.drawcoastlines();
    m.fillcontinents(color='grey');
    im1 = m.pcolormesh(gr.geolon,gr.geolat,df.zos[ind,:,:],shading='gouraud',
            vmin=-0.5,vmax=0.5,cmap='RdBu_r',latlon=True);
    cb = m.colorbar(im1,"left", size="5%", pad="2%",ticklocation='left')
    cb.ax.tick_params(labelsize=14)
    meridians = np.arange(-90.,-45.,10.)
    m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
    title=str(df.time[ind].coords)
    plt.title('SSH ' + title[-20:])
    ax = fig.add_subplot(122)
    m.drawcoastlines();
    m.fillcontinents(color='grey');
    im1 = m.pcolormesh(gr.geolon,gr.geolat,df.tos[ind,:,:],shading='gouraud',
            vmin=-1,vmax=28,cmap='copper_to_blue_r',latlon=True);
    cb = m.colorbar(im1,"right", size="5%", pad="2%")
    cb.ax.tick_params(labelsize=14)
    parallels = np.arange(25.,50.,5.)
    m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
    meridians = np.arange(-90.,-45.,10.)
    m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
    title=str(df.time[ind].coords)
    plt.title('SST ' + title[-20:])
    fname = 'gifs/GS_ssh_sst_' + str(itr).zfill(3) + '.png'
    plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)
    itr += 1
    plt.clf()
