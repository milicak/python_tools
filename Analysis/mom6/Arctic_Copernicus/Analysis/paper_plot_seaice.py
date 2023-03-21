import numpy as np

nsidc_year = np.arange(1992,2017)
nsidcarea = np.genfromtxt('/home/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis/NH_area_NSIDC_monthly.dat')
siarea_mnth = np.mean(nsidcarea,axis=0)
nsidcextent = np.genfromtxt('/home/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis/NH_extent_NSIDC_monthly.dat')
siextent_mnth = np.mean(nsidcextent,axis=0)

df1 = pd.read_excel('Sea_Ice_Index_Monthly_Data_by_Year_G02135_v3.xlsx', 'NH-Extent')
df1.rename( columns={'Unnamed: 0':'year'}, inplace=True )
df2 = pd.read_excel('Sea_Ice_Index_Monthly_Data_by_Year_G02135_v3.xlsx', 'NH-Area')
df2.rename( columns={'Unnamed: 0':'year'}, inplace=True )



dfc = xr.open_dataset('seaice_extent_area.nc')

# annual mean
siareaann = dfc.SI_area.groupby('Time.year').mean('Time')
siextentann = dfc.SI_extent.groupby('Time.year').mean('Time')
# monthly mean
siareamth = dfc.SI_area.groupby('Time.month').mean('Time')
siextentmth = dfce.SI_extent.groupby('Time.month').mean('Time')

# figure(1)
fig1, f1_axes = plt.subplots(figsize=(8,6),ncols=1, nrows=1,
                             constrained_layout=True, squeeze=False)
# f1_axes[0].plot(siareamth.month,siareamth.data*1e-12)
# f1_axes[0].plot(siareamth.month,siarea_mnth[:-1],'r')
# f1_axes[0].plot(siareamth.month,siextent_mnth[:-1],'k')
f1_axes[0,0].plot(df1['year'].iloc[18:40],df1['Annual'].iloc[18:40],label='NSIDC-extent')
f1_axes[0,0].plot(df2['year'].iloc[18:40],df2['Annual'].iloc[18:40],label='NSIDC-area')
f1_axes[0,0].plot(siextentann.year,siextentann.data*1e-12,label='MOM6-extent')
f1_axes[0,0].plot(siareaann.year,siareaann.data*1e-12,label='MOM6-area')
f1_axes[0,0].set_ylim(7,13);
f1_axes[0,0].set_xticklabels(["1990",'1995', "2000", "2005", "2010",'2015'], fontsize=14)
f1_axes[0,0].set_yticklabels(["7",'8', "9", "10", "11",'12','13'], fontsize=14)
f1_axes[0,0].set_ylabel(r'Sea-ice area [$10^6$ $\times$ km$^2$]',fontsize = 14.0)
f1_axes[0,0].set_xlabel(r'Year',fontsize = 14.0)
f1_axes[0,0].grid()
f1_axes[0,0].legend()

plt.savefig('paperfigs/seaice_obs_vs_model.png', bbox_inches='tight',format='png',dpi=300)





