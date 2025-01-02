from matplotlib import pyplot as plt, animation
gr = xr.open_dataset('ocean_geometry.nc')
df = xr.open_dataset('ocean_daily.nc')


ssh = np.zeros((df.time.shape[0],df.yh.shape[0],df.xh.shape[0]))

ssh[:,:,:] = np.copy(df.zos)


fig, ax = plt.subplots(figsize=(9,4))
cax = ax.pcolormesh(gr.geolon,gr.geolat,ssh[0,:,:], vmin=-.25, vmax=.25,
        shading='gouraud',cmap='BrBG_r')
plt.xlim(33.5,36.5);plt.ylim(35,37);
ax1 = plt.gca()
ax1.set_facecolor("grey")
fig.colorbar(cax)


def animate(i):
    cax.set_array(ssh[i,:,:].flatten())
    plt.xlim(33.5,36.5);plt.ylim(35,37);
    ax1 = plt.gca()
    ax1.set_facecolor("grey")
    title = str(df.time[i].coords)
    plt.title('Sea surface height ' + title[-5:] + ' minutes')


anim = animation.FuncAnimation(fig, animate, interval=140, frames=110)
anim.save('Antakya_Tsunami.gif')
# BEX4_tEJB3.U745
