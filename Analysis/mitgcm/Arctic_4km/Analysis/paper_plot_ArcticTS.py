import numpy as np

ds = xr.open_dataset('/archive/milicak/dataset/WOA13/woa13_decav_s00_04.nc',decode_times=False)
dt = xr.open_dataset('/archive/milicak/dataset/WOA13/woa13_decav_t00_04.nc',decode_times=False)

# Canada Basin
# lon = -140
# lat = 77
xind1 = 159
yind1 = 668

# Eurasia Basin
# lon = 38
# lat = 85
xind2 = 872
yind2 = 700

fig, axs = plt.subplots(1,2, sharey=True)
axs[0].plot(dt.t_an[0,:,yind1,xind1],-dt.depth,'b',label='Kanada Baseni')
axs[0].plot(dt.t_an[0,:,yind2,xind2],-dt.depth,'k',label='Avrasya Baseni')
axs[1].plot(ds.s_an[0,:,yind1,xind1],-dt.depth,'b',label='Kanada Baseni')
axs[1].plot(ds.s_an[0,:,yind2,xind2],-dt.depth,'k',label='Avrasya Baseni')

axs[1].yaxis.tick_right()
axs[0].tick_params(labelsize=14)
axs[1].tick_params(labelsize=14)
axs[0].spines['right'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[1].yaxis.set_ticks_position('right')
axs[0].set_ylabel('Derinlik [m]',fontsize=14)
axs[1].set_ylabel('Derinlik [m]',fontsize=14);
axs[1].yaxis.set_label_position("right")
axs[1].legend(loc='lower left')
axs[0].set_xlabel('Sıcaklık',fontsize=14)
axs[1].set_xlabel('Tuzluluk',fontsize=14)

plt.savefig('paperfigs/Arctic_real_TS_profiles.png', bbox_inches='tight',format='png',dpi=300)
