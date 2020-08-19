import numpy as np

# Cp*rho
cprho = 3996*1000

nsidc_year = np.arange(1992,2017)
nsidcarea = np.genfromtxt('NH_area_NSIDC_monthly.dat')
siarea_mnth = np.mean(nsidcarea,axis=0)
nsidcextent = np.genfromtxt('NH_extent_NSIDC_monthly.dat')
siextent_mnth = np.mean(nsidcextent,axis=0)

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/'

ctrl = root_folder + 'Exp02_0_seaice_area.nc'
ctrl_e = root_folder + 'Exp02_0_seaice_extent.nc'

dfc = xr.open_dataset(ctrl)
dfce = xr.open_dataset(ctrl_e)

# annual mean
siareaann = dfc.SIarea.groupby('time.year').mean('time')
siextentann = dfce.SI_extent.groupby('time.year').mean('time')
# monthly mean
siareamth = dfc.SIarea.groupby('time.month').mean('time')
siextentmth = dfce.SI_extent.groupby('time.month').mean('time')

# figure(1)
fig1, f1_axes = plt.subplots(figsize=(8,12),ncols=1, nrows=3, constrained_layout=True)
# f1_axes[0].plot(siareamth.month,siareamth.data*1e-12)
# f1_axes[0].plot(siareamth.month,siarea_mnth[:-1],'r')
# f1_axes[0].plot(siareamth.month,siextent_mnth[:-1],'k')
f1_axes[0].plot(siareaann.year,siareaann.data*1e-12,label='Ctrl')
f1_axes[0].plot(nsidc_year,nsidcarea[:,-1],'k',label='Obs')
f1_axes[0].set_ylim(8,12);
f1_axes[0].set_ylabel(r'Sea-ice area [$10^6$ $\times$ km$^2$]')
f1_axes[0].legend()

# Fram Strait Transport
ctrlFram = root_folder + 'Exp02_0_Fram_volume_transport.nc'
dfc_F = xr.open_dataset(ctrlFram)
timeref = dfc_F.time
dfc_F = dfc_F.volume_transport.groupby('time.year').mean('time')
# Davis Strait Transport
ctrlDavis = root_folder + 'Exp02_0_Davis_volume_transport.nc'
dfc_D = xr.open_dataset(ctrlDavis)
dfc_D = dfc_D.volume_transport.groupby('time.year').mean('time')
# Berin Strait Transport
ctrlBering = root_folder + 'Exp02_0_Bering_volume_transport.nc'
dfc_B = xr.open_dataset(ctrlBering)
dfc_B = dfc_B.volume_transport.groupby('time.year').mean('time')
# Barents Strait Transport
ctrlBarents = root_folder + 'Exp02_0_Barents_volume_transport.nc'
dfc_Br = xr.open_dataset(ctrlBarents)
dfc_Br['volume_transport'] = dfc_Br.v1_volume_transport + dfc_Br.v2_volume_transport - dfc_Br.u1_volume_transport
dfc_Br['time'] = timeref
dfc_Br = dfc_Br.volume_transport.groupby('time.year').mean('time')
# Into the Arctic
f1_axes[1].plot(dfc_F.year,dfc_F.data*1e-6,label='Fram')
f1_axes[1].plot(dfc_D.year,dfc_D.data*1e-6,label='Davis')
f1_axes[1].plot(dfc_B.year,-dfc_B.data*1e-6,label='Bering')
f1_axes[1].plot(dfc_Br.year,-dfc_Br.data*1e-6,label='Barents')
f1_axes[1].set_ylabel('Transport [Sv]')
f1_axes[1].legend()

# Fram Strait Heat Transport
ctrlFram = root_folder + 'Exp02_0_Fram_heat_transport.nc'
dfc_F = xr.open_dataset(ctrlFram)
dfc_F['heat_transport'] = 0.6*dfc_F.heat_transport
dfc_F = dfc_F.heat_transport.groupby('time.year').mean('time')
# Davis Strait  Heat Transport
ctrlDavis = root_folder + 'Exp02_0_Davis_heat_transport.nc'
dfc_D = xr.open_dataset(ctrlDavis)
dfc_D = dfc_D.heat_transport.groupby('time.year').mean('time')
# Berin Strait Heat Transport
ctrlBering = root_folder + 'Exp02_0_Bering_heat_transport.nc'
dfc_B = xr.open_dataset(ctrlBering)
dfc_B = dfc_B.heat_transport.groupby('time.year').mean('time')
# Barents Strait Heat Transport
ctrlBarents = root_folder + 'Exp02_0_Barents_heat_transport.nc'
dfc_Br = xr.open_dataset(ctrlBarents)
dfc_Br['heat_transport'] = dfc_Br.v1_heat_transport + dfc_Br.v2_heat_transport - dfc_Br.u1_heat_transport
dfc_Br['time'] = timeref
dfc_Br = dfc_Br.heat_transport.groupby('time.year').mean('time')
# Into the Arctic
f1_axes[2].plot(dfc_F.year,dfc_F.data*cprho*1e-12,label='Fram')
f1_axes[2].plot(dfc_D.year,dfc_D.data*cprho*1e-12,label='Davis')
f1_axes[2].plot(dfc_B.year,-dfc_B.data*cprho*1e-12,label='Bering')
f1_axes[2].plot(dfc_Br.year,-dfc_Br.data*cprho*1e-12,label='Barents')
f1_axes[2].set_ylabel('Heat Transport [TW]')
f1_axes[2].set_xlabel('Year')
f1_axes[2].set_ylim(-5,130);
# f1_axes[2].legend()
f1_axes[2].legend(loc='upper right')

plt.savefig('paperfigs/ctrl_seaice_volume_heat.png', bbox_inches='tight',format='png',dpi=300)





