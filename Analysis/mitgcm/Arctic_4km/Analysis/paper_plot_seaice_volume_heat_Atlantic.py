import numpy as np

# Cp*rho
cprho = 3996*1000

nsidc_year = np.arange(1992,2017)
nsidcarea = np.genfromtxt('NH_area_NSIDC_monthly.dat')
siarea_mnth = np.mean(nsidcarea,axis=0)
nsidcextent = np.genfromtxt('NH_extent_NSIDC_monthly.dat')
siextent_mnth = np.mean(nsidcextent,axis=0)

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/'

# ctrl = root_folder + 'Exp02_0_seaice_area.nc'
ctrl = root_folder + 'Exp02_0_seaice_extent.nc'
dfc = xr.open_dataset(ctrl)
# warm1 = root_folder + 'Exp02_1_seaice_area.nc'
warm1 = root_folder + 'Exp02_1_seaice_extent.nc'
dfw1 = xr.open_dataset(warm1)

# annual mean
# siareaann_w1 = dfw1.SI_area.groupby('time.year').mean('time')-dfc.SIarea.groupby('time.year').mean('time')
siextentann_w1 = dfw1.SI_extent.groupby('time.year').mean('time')-dfc.SI_extent.groupby('time.year').mean('time')

# control volume transport
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

# Atlantic Volume Transport
W1Fram = root_folder + 'Exp02_1_Fram_volume_transport.nc'
dfw_F = xr.open_dataset(W1Fram)
timeref = dfw_F.time
dfw_F = dfw_F.volume_transport.groupby('time.year').mean('time')
# Davis Strait Transport
W1Davis = root_folder + 'Exp02_1_Davis_volume_transport.nc'
dfw_D = xr.open_dataset(W1Davis)
dfw_D = dfw_D.volume_transport.groupby('time.year').mean('time')
# Berin Strait Transport
W1Bering = root_folder + 'Exp02_1_Bering_volume_transport.nc'
dfw_B = xr.open_dataset(W1Bering)
dfw_B = dfw_B.volume_transport.groupby('time.year').mean('time')
# Barents Strait Transport
W1Barents = root_folder + 'Exp02_1_Barents_volume_transport.nc'
dfw_Br = xr.open_dataset(W1Barents)
dfw_Br['volume_transport'] = dfw_Br.v1_volume_transport + dfw_Br.v2_volume_transport - dfw_Br.u1_volume_transport
dfw_Br['time'] = timeref
dfw_Br = dfw_Br.volume_transport.groupby('time.year').mean('time')

# Fram Strait Heat Transport
ctrlhFram = root_folder + 'Exp02_0_Fram_heat_transport.nc'
dfch_F = xr.open_dataset(ctrlhFram)
dfch_F = dfch_F.heat_transport.groupby('time.year').mean('time')
# Davis Strait  Heat Transport
ctrlhDavis = root_folder + 'Exp02_0_Davis_heat_transport.nc'
dfch_D = xr.open_dataset(ctrlhDavis)
dfch_D = dfch_D.heat_transport.groupby('time.year').mean('time')
# Berin Strait Heat Transport
ctrlhBering = root_folder + 'Exp02_0_Bering_heat_transport.nc'
dfch_B = xr.open_dataset(ctrlhBering)
dfch_B = dfch_B.heat_transport.groupby('time.year').mean('time')
# Barents Strait Heat Transport
ctrlhBarents = root_folder + 'Exp02_0_Barents_heat_transport.nc'
dfch_Br = xr.open_dataset(ctrlhBarents)
dfch_Br['heat_transport'] = dfch_Br.v1_heat_transport + dfch_Br.v2_heat_transport - dfch_Br.u1_heat_transport
dfch_Br['time'] = timeref
dfch_Br = dfch_Br.heat_transport.groupby('time.year').mean('time')

# Atlantic Simulation Heat Transport
W1hFram = root_folder + 'Exp02_1_Fram_heat_transport.nc'
dfwh_F = xr.open_dataset(W1hFram)
dfwh_F = dfwh_F.heat_transport.groupby('time.year').mean('time')
# Davis Strait  Heat Transport
W1hDavis = root_folder + 'Exp02_1_Davis_heat_transport.nc'
dfwh_D = xr.open_dataset(W1hDavis)
dfwh_D = dfwh_D.heat_transport.groupby('time.year').mean('time')
# Berin Strait Heat Transport
W1hBering = root_folder + 'Exp02_1_Bering_heat_transport.nc'
dfwh_B = xr.open_dataset(W1hBering)
dfwh_B = dfwh_B.heat_transport.groupby('time.year').mean('time')
# Barents Strait Heat Transport
W1hBarents = root_folder + 'Exp02_1_Barents_heat_transport.nc'
dfwh_Br = xr.open_dataset(W1hBarents)
dfwh_Br['heat_transport'] = dfwh_Br.v1_heat_transport + dfwh_Br.v2_heat_transport - dfwh_Br.u1_heat_transport
dfwh_Br['time'] = timeref
dfwh_Br = dfwh_Br.heat_transport.groupby('time.year').mean('time')

# figure(1)
fig1, f1_axes = plt.subplots(figsize=(8,12),ncols=1, nrows=3, constrained_layout=True)
f1_axes[0].plot(siareaann.year,siareaann_w1.data*1e-12,label='Ctrl')
# f1_axes[0].set_ylim(8,12);
f1_axes[0].set_ylabel(r'Sea-ice area [$10^6$ $\times$ km$^2$]')
f1_axes[0].legend()

# Fram Strait Transport
# Into the Arctic
f1_axes[1].plot(dfc_F.year,(dfw_F.data-dfc_F.data)*1e-6,label='Fram')
f1_axes[1].plot(dfc_Br.year,-(dfw_Br.data-dfc_Br.data)*1e-6,label='Barents')
f1_axes[1].set_ylabel('Transport [Sv]')
f1_axes[1].legend()

# Into the Arctic
f1_axes[2].plot(dfc_F.year,(dfwh_F.data-dfch_F.data)*cprho*1e-12,label='Fram')
f1_axes[2].plot(dfc_Br.year,-(dfwh_Br.data-dfch_Br.data)*cprho*1e-12,label='Barents')
f1_axes[2].set_ylabel('Heat Transport [TW]')
f1_axes[2].set_xlabel('Year')
f1_axes[2].set_ylim(-5,130);
# f1_axes[2].legend()
f1_axes[2].legend(loc='lower right')

plt.savefig('paperfigs/Atlantic_warm_seaice_volume_heat.png', bbox_inches='tight',format='png',dpi=300)





