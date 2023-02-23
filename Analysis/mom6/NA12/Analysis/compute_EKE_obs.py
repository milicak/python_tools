import numpy as np

ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/SSHUGVG/global/*.nc'))

df = xr.open_mfdataset(ls1)
# The geostrophic currents are derived
# from sla (geostrophic velocities anomalies, ugosa and vgosa variables) and
# from adt (absolute geostrophic velicities, ugos and vgos variables).

umean = df.ugos.mean('time')
vmean = df.vgos.mean('time')
upr = df.ugos - umean
vpr = df.vgos - vmean

EKE = upr**2+vpr**2
EKE = EKE.mean('time')


ds = EKE.to_dataset(name='EKE')
ds.to_netcdf('/Volumes/A1/workdir/milicak/datasets/SSHUGVG/global_EKE.nc')
