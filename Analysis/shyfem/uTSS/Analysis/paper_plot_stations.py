df=pd.read_csv('obs_ctd/2017_07_K0_obs.csv')

dfm=pd.read_csv('obs_ctd/2017_07_K0_uTSS.csv')



fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(dfm.Salt_uTSS,dfm.zr_uTSS);
ax1.plot(df.Salt_obs,df.zr_obs);
ax2.plot(dfm.Temp_uTSS,dfm.zr_uTSS);
ax2.plot(df.Temp_obs,df.zr_obs);

ax1.set_ylabel('Depth [m]');
ax1.set_xlabel('Salinity [psu]');
ax2.set_xlabel('Temperature [$^\circ$C]');
ax2.legend(('model','obs'));

plt.savefig('paperfigs/K0_07_2017_model_vs_obs.png', bbox_inches='tight',format='png', dpi=300)

