plt.figure(figsize=(8,4))
#plt.clf();
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#gl.xlabels_top = False
#     ...: gl.ylabels_left = False
#     ...: gl.xlines = False
#     ...: gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
#     ...: gl.xformatter = LONGITUDE_FORMATTER
#     ...: gl.yformatter = LATITUDE_FORMATTER
#     ...: gl.xlabel_style = {'size': 15, 'color': 'gray'}
#     ...: gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

plt.pcolormesh(data.lon,data.lat,dnm-dnmref,vmin=-3,vmax=3,cmap='RdYlBu_r',transform=ccrs.PlateCarree())
plt.colorbar(ax=ax, shrink=.75)
plt.tight_layout()
