import numpy
import xarray as xr


df=xr.open_dataset('temp_salt_thalweg_sections.nc')

fig, axes = plt.subplots(figsize=(14,6))                   
ax1 = plt.subplot2grid(shape=(1,2),loc=(0,0), colspan=1)
ax2 = plt.subplot2grid(shape=(1,2),loc=(0,1), colspan=1)

plt.tight_layout()                                         

c1 = ax1.pcolormesh(df.dist,-df.zr,np.transpose(ma.masked_invalid(df.salt_thalweg))
              ,cmap='nasa_rainbow',vmin=17.5,vmax=38.5);
ax1.set_ylim(-500,0);
ax1.invert_xaxis()
ax1.set_xticklabels(['400', '350', '300', '250', '200', '150', '100', '50',
                     '0'], fontsize=14);
axpos = ax1.get_position()                                           
cbar_ax = fig.add_axes([axpos.x1-0.4,axpos.y0+0.05,0.03,axpos.height*0.6]) 
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='right')           
cbar.set_label('psu',rotation=0,y=1.07,labelpad=-45);
ax1.set_yticklabels(ax1.get_yticks(), rotation=0, fontsize=14);

c2 = ax2.pcolormesh(df.dist,-df.zr,np.transpose(ma.masked_invalid(df.temp_thalweg))
              ,cmap='needJet2',vmin=9,vmax=18);
ax2.set_ylim(-500,0);
ax2.invert_xaxis()
ax2.set_xticklabels(['400', '350', '300', '250', '200', '150', '100', '50',
                     '0'], fontsize=14);
axpos = ax2.get_position()                                           
cbar_ax = fig.add_axes([axpos.x1-0.4,axpos.y0+0.05,0.03,axpos.height*0.6]) 
cbar = fig.colorbar(c2, cax=cbar_ax, ticklocation='right')           
cbar.set_label('$^\circ$C',rotation=0,y=1.07,labelpad=-45);
ax2.set_yticklabels(ax2.get_yticks(), rotation=0, fontsize=14);

ax1.set_xlabel('Distance (km)', fontsize=14);
ax1.set_ylabel('Depth (m)', fontsize=14);
ax2.set_xlabel('Distance (km)', fontsize=14);
ax2.yaxis.tick_right() 

from matplotlib.ticker import FuncFormatter 
def minus_formatter(x, pos):
    return unicode(x).replace('-', u'\u22122')

def math_formatter(x, pos):
    return "$%s$" % x

ax1.yaxis.set_major_formatter(FuncFormatter(math_formatter)) 
ax2.yaxis.set_major_formatter(FuncFormatter(math_formatter)) 

# plt.savefig('paperfigs/Salt_Temp_mean_thalweg.png', bbox_inches='tight',format='png',dpi=300)
