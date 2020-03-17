import gcsfs
from matplotlib import path

coeff = 1.0
# cesm
dfct = xr.open_dataset('CESM_tempSP_ssp585.nc');
dfcs = xr.open_dataset('CESM_saltSP_ssp585.nc');
dfca = xr.open_dataset('CESM_amoc_ssp585.nc');

plt.figure()
plt.plot(-0.2*(dfct.sstSP-dfct.sstSP[0])+0.8*(dfcs.sssSP-dfcs.sssSP[0]),dfca.msftmz-coeff*dfca.msftmz[0],'g*')

# access
dfat = xr.open_dataset('ACCESS_tempSP_ssp585.nc');
dfas = xr.open_dataset('ACCESS_saltSP_ssp585.nc');
dfaa = xr.open_dataset('ACCESS_amoc_ssp585.nc');
# plt.figure()
plt.plot(-0.2*(dfat.sstSP-dfat.sstSP[0])+0.8*(dfas.sssSP-dfas.sssSP[0]),dfaa.msftmz-coeff*dfaa.msftmz[0],'k*')

# mpi
dfmt = xr.open_dataset('MPI_tempSP_ssp585.nc');
dfms = xr.open_dataset('MPI_saltSP_ssp585.nc');
dfma = xr.open_dataset('MPI_amoc_ssp585.nc');
# plt.figure()
plt.plot(-0.2*(dfmt.sstSP-dfmt.sstSP[0])+0.8*(dfms.sssSP-dfms.sssSP[0]),dfma.msftmz-coeff*dfma.msftmz[0],'r*')

# noresm
dfnt = xr.open_dataset('NorESM_MM_tempSP_ssp585.nc');
dfns = xr.open_dataset('NorESM_MM_saltSP_ssp585.nc');
dfna = xr.open_dataset('NorESM_MM_amoc_ssp585.nc');
# plt.figure()
plt.plot(-0.2*(dfnt.sstSP-dfnt.sstSP[0])+0.8*(dfns.sssSP-dfns.sssSP[0]),dfna.msftmz-coeff*dfna.msftmz[0],'b*')

dfKt = xr.open_dataset('CanESM_tempSP_ssp585.nc');
dfKs = xr.open_dataset('CanESM_saltSP_ssp585.nc');
dfKa = xr.open_dataset('CanESM_amoc_ssp585.nc');
dfKa = dfKa.sel(year=slice(2015,2100))
# plt.figure()
plt.plot(-0.2*(dfKt.sstSP-dfKt.sstSP[0])+0.8*(dfKs.sssSP-dfKs.sssSP[0]),dfKa.msftmz-coeff*dfKa.msftmz[0],'m*')



dfntt = xr.open_dataset('NorESM_MM_tempST_ssp585.nc');
dfnst = xr.open_dataset('NorESM_MM_saltST_ssp585.nc');
dfmtt = xr.open_dataset('MPI_tempST_ssp585.nc');
dfmst = xr.open_dataset('MPI_saltST_ssp585.nc');
dfatt = xr.open_dataset('ACCESS_tempST_ssp585.nc');
dfast = xr.open_dataset('ACCESS_saltST_ssp585.nc');
dfctt = xr.open_dataset('CESM_tempST_ssp585.nc');
dfcst = xr.open_dataset('CESM_saltST_ssp585.nc');
dfKtt = xr.open_dataset('CanESM_tempST_ssp585.nc');
dfKst = xr.open_dataset('CanESM_saltST_ssp585.nc');
