import numpy as np
from scipy.io import loadmat

# var = loadmat('/archive/milicak/MITgcm_c65/Projects/thermobaricity/Exp01.0/Exp01.0_31680_eps_temp.mat')
# var = loadmat('/archive/milicak/MITgcm_c65/Projects/thermobaricity/Exp01.1/Exp01.1_52560_eps_temp.mat')
var = loadmat('/archive/milicak/MITgcm_c65/Projects/thermobaricity/Exp01.2/Exp01.2_288000_eps_temp.mat')
# var = loadmat('/archive/milicak/MITgcm_c65/Projects/thermobaricity/Exp01.3/Exp01.3_86400_eps_temp.mat')
# savename = 'paperfigs/snapshots_Exp01_0.png'
# savename = 'paperfigs/snapshots_Exp01_1.png'
savename = 'paperfigs/snapshots_Exp01_2.png'
# savename = 'paperfigs/snapshots_Exp01_3.png'

x = var['x']
z = var['Z']
epsilon = var['epsilon_rho']
temp = var['temp']

epsilon[:,:,2300::] = 1e-10;

fig, axs = plt.subplots(1,2,figsize=(15,6))
t1_plot = axs[0].pcolormesh(x[:,0],-z[:-1,:],np.transpose(temp[:,0,:]),cmap='nice_gfdl',rasterized=True);
cb = plt.colorbar(t1_plot,ax=axs[0]);
cb.ax.tick_params(labelsize=14)
t1_plot.set_clim(-2,1.5)
axs[0].set_ylabel('Depth [m]', fontsize=14)
axs[0].set_yticklabels([-3000,-2500,-2000,-1500,-1000,-500,0], fontsize=14)
axs[0].set_xticklabels([0,500,1000,1500,2000,2500,3000,3500,4000], fontsize=14)
axs[0].set_xlabel('x [m]', fontsize=14);

e1_plot = axs[1].pcolormesh(x[:,0],-z[:-1,:],np.transpose(np.log10(epsilon[:,1,:])),cmap='needJet2',rasterized=True);
cb = plt.colorbar(e1_plot,ax=axs[1]);
cb.ax.tick_params(labelsize=14)
e1_plot.set_clim(-7,-3)
axs[1].set_ylabel('Depth [m]', fontsize=14)
axs[1].set_yticklabels([-3000,-2500,-2000,-1500,-1000,-500,0], fontsize=14)
axs[1].set_xticklabels([0,500,1000,1500,2000,2500,3000,3500,4000], fontsize=14)
axs[1].set_xlabel('x [m]', fontsize=14);

fig.canvas.draw()
plt.tight_layout()

plt.savefig(savename, bbox_inches='tight',format='png',dpi=300)
