import numpy as np
import xarray as xr



prenames = ['2017_08_', '2018_01_', '2018_04_', '2018_08_']


fnames = ['obs_ctd_netcdf_athena/2017_08_',
          'obs_ctd_netcdf/2018_01_',
          'obs_ctd_netcdf/2018_04_',
          'obs_ctd_netcdf_athena/2018_08_']

dfb1 = pd.read_csv('obs_ctd_netcdf_athena/2017_08_bias.csv')
dfb2 = pd.read_csv('obs_ctd_netcdf/2018_01_bias.csv')
dfb3 = pd.read_csv('obs_ctd_netcdf/2018_04_bias.csv')
dfb4 = pd.read_csv('obs_ctd_netcdf_athena/2018_08_bias.csv')

dfr1 = pd.read_csv('obs_ctd_netcdf_athena/2017_08_rms.csv')
dfr2 = pd.read_csv('obs_ctd_netcdf/2018_01_rms.csv')
dfr3 = pd.read_csv('obs_ctd_netcdf_athena/2018_04_rms.csv')
dfr4 = pd.read_csv('obs_ctd_netcdf_athena/2018_08_rms.csv')

def plot_two_data(df1,df2,df3,df4,var1='Salt',var2='Temp'):
    fig, axs = plt.subplots(1,2, sharey=True)
    axs[0].plot(df1[var1]*0.8,df1.zr_obs,'b',label='2017-08')
    axs[0].plot(df2[var1]*0.8,df2.zr_obs,'r',label='2018-01')
    axs[0].plot(df3[var1]*0.8,df3.zr_obs,'g',label='2018-04')
    axs[0].plot(df4[var1]*0.8,df4.zr_obs,'k',label='2018-08')
    axs[1].plot(df1[var2]*0.8,df1.zr_obs,'b',label='2017-08')
    axs[1].plot(df2[var2]*0.8,df2.zr_obs,'r',label='2018-01')
    axs[1].plot(df3[var2]*0.8,df3.zr_obs,'g',label='2018-04')
    axs[1].plot(df4[var2]*0.8,df4.zr_obs,'k',label='2018-08')
    axs[1].yaxis.tick_right()
    axs[0].tick_params(labelsize=14)
    axs[1].tick_params(labelsize=14)
    axs[0].spines['right'].set_visible(False)
    axs[1].spines['left'].set_visible(False)
    axs[1].spines['top'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    axs[1].yaxis.set_ticks_position('right')
    axs[0].set_ylabel('Depth (m)',fontsize=14)
    axs[1].set_ylabel('Depth (m)',fontsize=14);
    axs[1].yaxis.set_label_position("right")

    return axs


axs = plot_two_data(dfb1,dfb2,dfb3,dfb4,'Salt_bias','Temp_bias')
axs[0].set_xlabel('Salinity bias (psu)',fontsize=14)
axs[1].set_xlabel('Temperature bias ($^\circ$C)',fontsize=14)
axs[0].set_xlim(-2,5)
axs[1].set_xlim(-2.5,2.5)
axs[0].legend(['2017-08','2018-01','2018-04','2018-08'],fontsize=14,frameon=False)
plt.savefig('paperfigs/salinity_temperature_bias_all.png', bbox_inches='tight',format='png',dpi=300)

axs = plot_two_data(dfr1,dfr2,dfr3,dfr4,'Salt_rms','Temp_rms')
axs[0].set_xlabel('Salinity rms (psu)',fontsize=14)
axs[1].set_xlabel('Temperature rms ($^\circ$C)',fontsize=14)
axs[0].set_xlim(0,5)
axs[1].set_xlim(0,5)
axs[0].legend(['2017-08','2018-01','2018-04','2018-08'],fontsize=14,frameon=False)
plt.savefig('paperfigs/salinity_temperature_rms_all.png', bbox_inches='tight',format='png',dpi=300)

