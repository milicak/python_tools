import numpy as np


def diff_norm(df1,df2,var='SI_extent'):
    df = 100 * (df1[var].groupby('time.year').mean('time')
          - df2[var].groupby('time.year').mean('time'))/(df2[var].groupby('time.year').mean('time'))
    return df

def diff_dim(df1,df2,var='SI_extent'):
    df = (df1[var].groupby('time.year').mean('time')
          - df2[var].groupby('time.year').mean('time'))
    return df

# Cp*rho
cprho = 3996*1025

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/'

ctrl = root_folder + 'Exp02_0_seaice_extent.nc'
dfc = xr.open_dataset(ctrl)
warm1 = root_folder + 'Exp02_1_seaice_extent.nc'
dfw1 = xr.open_dataset(warm1)
warm3 = root_folder + 'Exp02_3_seaice_extent.nc'
dfw3 = xr.open_dataset(warm3)

ctrl = root_folder + 'Exp02_0_seaice_extent_Barents.nc'
dfcB = xr.open_dataset(ctrl)
warm1 = root_folder + 'Exp02_1_seaice_extent_Barents.nc'
dfw1B = xr.open_dataset(warm1)
warm3 = root_folder + 'Exp02_3_seaice_extent_Barents.nc'
dfw3B = xr.open_dataset(warm3)

# annual mean
SI_ext_NH_w1 = diff_norm(dfw1,dfc,'SI_extent')
SI_ext_NH_w3 = diff_norm(dfw3,dfc,'SI_extent')
SI_ext_Br_w1 = diff_norm(dfw1B,dfcB,'SI_extent')
SI_ext_Br_w3 = diff_norm(dfw3B,dfcB,'SI_extent')

# control volume transport
# Fram Strait Transport
ctrlFram = root_folder + 'Exp02_0_Fram_volume_transport.nc'
dfc_F = xr.open_dataset(ctrlFram)
timeref = dfc_F.time
warm1Fram = root_folder + 'Exp02_1_Fram_volume_transport.nc'
dfw1_F = xr.open_dataset(warm1Fram)
warm3Fram = root_folder + 'Exp02_3_Fram_volume_transport.nc'
dfw3_F = xr.open_dataset(warm3Fram)
dfw3_F['time'] = timeref
Vol_tr_F_w1 = diff_dim(dfw1_F,dfc_F,'volume_transport')
Vol_tr_F_w3 = diff_dim(dfw3_F,dfc_F,'volume_transport')
# Davis Strait Transport
ctrlDavis = root_folder + 'Exp02_0_Davis_volume_transport.nc'
dfc_D = xr.open_dataset(ctrlDavis)
warm1Davis = root_folder + 'Exp02_1_Davis_volume_transport.nc'
dfw1_D = xr.open_dataset(warm1Davis)
warm3Davis = root_folder + 'Exp02_3_Davis_volume_transport.nc'
dfw3_D = xr.open_dataset(warm3Davis)
dfw3_D['time'] = timeref
Vol_tr_D_w1 = diff_dim(dfw1_D,dfc_D,'volume_transport')
Vol_tr_D_w3 = diff_dim(dfw3_D,dfc_D,'volume_transport')
# Bering Strait Transport
ctrlBering = root_folder + 'Exp02_0_Bering_volume_transport.nc'
dfc_B = xr.open_dataset(ctrlBering)
warm1Bering = root_folder + 'Exp02_1_Bering_volume_transport.nc'
dfw1_B = xr.open_dataset(warm1Bering)
warm3Bering = root_folder + 'Exp02_3_Bering_volume_transport.nc'
dfw3_B = xr.open_dataset(warm3Bering)
dfw3_B['time'] = timeref
Vol_tr_B_w1 = diff_dim(dfw1_B,dfc_B,'volume_transport')
Vol_tr_B_w3 = diff_dim(dfw3_B,dfc_B,'volume_transport')
# Barents Strait Transport
ctrlBarents = root_folder + 'Exp02_0_Barents_volume_transport.nc'
dfc_Br = xr.open_dataset(ctrlBarents)
dfc_Br['volume_transport'] = dfc_Br.v1_volume_transport + dfc_Br.v2_volume_transport - dfc_Br.u1_volume_transport
dfc_Br['time'] = timeref
warm1Barents = root_folder + 'Exp02_1_Barents_volume_transport.nc'
dfw1_Br = xr.open_dataset(warm1Barents)
dfw1_Br['volume_transport'] = dfw1_Br.v1_volume_transport + dfw1_Br.v2_volume_transport - dfw1_Br.u1_volume_transport
dfw1_Br['time'] = timeref
warm3Barents = root_folder + 'Exp02_3_Barents_volume_transport.nc'
dfw3_Br = xr.open_dataset(warm3Barents)
dfw3_Br['volume_transport'] = dfw3_Br.v1_volume_transport + dfw3_Br.v2_volume_transport - dfw3_Br.u1_volume_transport
dfw3_Br['time'] = timeref
Vol_tr_Br_w1 = diff_dim(dfw1_Br,dfc_Br,'volume_transport')
Vol_tr_Br_w3 = diff_dim(dfw3_Br,dfc_Br,'volume_transport')


# control heat transport
# Fram Strait Transport
ctrlFram = root_folder + 'Exp02_0_Fram_heat_transport.nc'
dfc_F = xr.open_dataset(ctrlFram)
timeref = dfc_F.time
warm1Fram = root_folder + 'Exp02_1_Fram_heat_transport.nc'
dfw1_F = xr.open_dataset(warm1Fram)
warm3Fram = root_folder + 'Exp02_3_Fram_heat_transport.nc'
dfw3_F = xr.open_dataset(warm3Fram)
dfw3_F['time'] = timeref
HT_tr_F_w1 = diff_dim(dfw1_F,dfc_F,'heat_transport')
HT_tr_F_w3 = diff_dim(dfw3_F,dfc_F,'heat_transport')
# Davis Strait Transport
ctrlDavis = root_folder + 'Exp02_0_Davis_heat_transport.nc'
dfc_D = xr.open_dataset(ctrlDavis)
warm1Davis = root_folder + 'Exp02_1_Davis_heat_transport.nc'
dfw1_D = xr.open_dataset(warm1Davis)
warm3Davis = root_folder + 'Exp02_3_Davis_heat_transport.nc'
dfw3_D = xr.open_dataset(warm3Davis)
dfw3_D['time'] = timeref
HT_tr_D_w1 = diff_dim(dfw1_D,dfc_D,'heat_transport')
HT_tr_D_w3 = diff_dim(dfw3_D,dfc_D,'heat_transport')
# Bering Strait Transport
ctrlBering = root_folder + 'Exp02_0_Bering_heat_transport.nc'
dfc_B = xr.open_dataset(ctrlBering)
warm1Bering = root_folder + 'Exp02_1_Bering_heat_transport.nc'
dfw1_B = xr.open_dataset(warm1Bering)
warm3Bering = root_folder + 'Exp02_3_Bering_heat_transport.nc'
dfw3_B = xr.open_dataset(warm3Bering)
dfw3_B['time'] = timeref
HT_tr_B_w1 = diff_dim(dfw1_B,dfc_B,'heat_transport')
HT_tr_B_w3 = diff_dim(dfw3_B,dfc_B,'heat_transport')
# Barents Strait Transport
ctrlBarents = root_folder + 'Exp02_0_Barents_heat_transport.nc'
dfc_Br = xr.open_dataset(ctrlBarents)
dfc_Br['heat_transport'] = dfc_Br.v1_heat_transport + dfc_Br.v2_heat_transport - dfc_Br.u1_heat_transport
dfc_Br['time'] = timeref
warm1Barents = root_folder + 'Exp02_1_Barents_heat_transport.nc'
dfw1_Br = xr.open_dataset(warm1Barents)
dfw1_Br['heat_transport'] = dfw1_Br.v1_heat_transport + dfw1_Br.v2_heat_transport - dfw1_Br.u1_heat_transport
dfw1_Br['time'] = timeref
warm3Barents = root_folder + 'Exp02_3_Barents_heat_transport.nc'
dfw3_Br = xr.open_dataset(warm3Barents)
dfw3_Br['heat_transport'] = dfw3_Br.v1_heat_transport + dfw3_Br.v2_heat_transport - dfw3_Br.u1_heat_transport
dfw3_Br['time'] = timeref
HT_tr_Br_w1 = diff_dim(dfw1_Br,dfc_Br,'heat_transport')
HT_tr_Br_w3 = diff_dim(dfw3_Br,dfc_Br,'heat_transport')

linecolors = (plt.rcParams['axes.prop_cycle'].by_key()['color'])

fig1, f1_axes = plt.subplots(figsize=(8,12),ncols=1, nrows=3, constrained_layout=True)
f1_axes[0].plot(SI_ext_NH_w1.year,SI_ext_NH_w1.data,
                color=linecolors[0],label='NH-Atlantic1')
f1_axes[0].plot(SI_ext_NH_w3.year,SI_ext_NH_w3.data,
                color=linecolors[1],label='NH-Atlantic2')
f1_axes[0].plot(SI_ext_NH_w1.year,SI_ext_Br_w1.data,
                color=linecolors[0],linestyle='dashed',label='Ba-Atl1')
f1_axes[0].plot(SI_ext_NH_w3.year,SI_ext_Br_w3.data,
                color=linecolors[1],linestyle='dashed',label='Ba-Atl2')
f1_axes[0].set_ylim(-8,8);
f1_axes[0].set_ylabel(r'Change in Sea-ice area [$\%$]',fontsize=16)
f1_axes[0].legend(fontsize=14,loc='upper right',frameon=False,ncol=2)
f1_axes[0].tick_params(labelsize=16)
# f1_axes[0].grid()
f1_axes[0].text(1991.1,-6.8,'a)',fontsize=16)

# Into the Arctic
f1_axes[1].plot(Vol_tr_F_w1.year,Vol_tr_F_w1.data*1e-6,color=linecolors[0],label='Fram')
f1_axes[1].plot(Vol_tr_F_w3.year,Vol_tr_F_w3.data*1e-6,color=linecolors[0],linestyle='dashed',label='Fram')
f1_axes[1].plot(Vol_tr_D_w1.year,Vol_tr_D_w1.data*1e-6,color=linecolors[1],label='Davis')
f1_axes[1].plot(Vol_tr_D_w3.year,Vol_tr_D_w3.data*1e-6,color=linecolors[1],linestyle='dashed',label='Davis')
f1_axes[1].plot(Vol_tr_B_w1.year,-Vol_tr_B_w1.data*1e-6,color=linecolors[2],label='Bering')
f1_axes[1].plot(Vol_tr_B_w3.year,-Vol_tr_B_w3.data*1e-6,color=linecolors[2],linestyle='dashed',label='Bering')
f1_axes[1].plot(Vol_tr_Br_w1.year,-Vol_tr_Br_w1.data*1e-6,color=linecolors[3],label='Barents')
f1_axes[1].plot(Vol_tr_Br_w3.year,-Vol_tr_Br_w3.data*1e-6,color=linecolors[3],linestyle='dashed',label='Barents')
f1_axes[1].set_ylabel('Transport [Sv]',fontsize=16)
f1_axes[1].tick_params(labelsize=16)
# f1_axes[1].grid()
f1_axes[1].text(1991.1,-0.6,'b)',fontsize=16)


f1_axes[2].plot(HT_tr_F_w1.year,cprho*HT_tr_F_w1.data*1e-12,color=linecolors[0],label='Fr-Atl1')
f1_axes[2].plot(HT_tr_F_w3.year,cprho*HT_tr_F_w3.data*1e-12,color=linecolors[0],linestyle='dashed',label='Fr-Atl2')
f1_axes[2].plot(HT_tr_D_w1.year,cprho*HT_tr_D_w1.data*1e-12,color=linecolors[1],label='Da-Atl1')
f1_axes[2].plot(HT_tr_D_w3.year,cprho*HT_tr_D_w3.data*1e-12,color=linecolors[1],linestyle='dashed',label='Da-Atl2')
f1_axes[2].plot(HT_tr_B_w1.year,-cprho*HT_tr_B_w1.data*1e-12,color=linecolors[2],label='Be-Atl1')
f1_axes[2].plot(HT_tr_B_w3.year,-cprho*HT_tr_B_w3.data*1e-12,color=linecolors[2],linestyle='dashed',label='Be-Atl2')
f1_axes[2].plot(HT_tr_Br_w1.year,-cprho*HT_tr_Br_w1.data*1e-12,color=linecolors[3],label='Ba-Atl1')
f1_axes[2].plot(HT_tr_Br_w3.year,-cprho*HT_tr_Br_w3.data*1e-12,color=linecolors[3],linestyle='dashed',label='Ba-Atl2')

f1_axes[2].set_ylabel('Heat Transport [TW]',fontsize=16)
# f1_axes[2].grid()
f1_axes[2].set_xlabel('Year',fontsize=16)
f1_axes[2].set_ylim(-40,40);
f1_axes[2].legend(fontsize=14,loc='lower center',frameon=False,ncol=4)
f1_axes[2].tick_params(labelsize=16)
f1_axes[2].text(1991.1,-38,'c)',fontsize=16);

plt.savefig('paperfigs/Atlantic_seaice_volume_heat.png', bbox_inches='tight',format='png',dpi=300)

