import numpy as np

if __name__ == '__main__':
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=4, memory_limit=25e9, local_dir='/cluster/work/users/cgu025/tmp_annual_mean/')
    client = Client(cluster)


root_folder = '/tos-project1/NS9252K/CMIP6/'
omipn = 'omip1'
model = 'GFDL-CM4'
scn = 'r1i1p1f1'

fnames = root_folder+omipn+'/'+ model+'/'+scn+'/thetao_Omon_'+model+'_'+omipn+'_'+scn+'_gn*'
# fnames = root_folder+omipn+'/'+ model+'/'+scn+'/so_Omon_'+model+'_'+omipn+'_'+scn+'_gn*'

list=sorted(glob.glob(fnames))

CHUNKS = {'lev': 1}
# CHUNKS = {'y': 270, 'x': 360}


# df=xr.open_mfdataset(list,chunks=CHUNKS)
df=xr.open_mfdataset(list)
# GFDL indices
# Fram Strait
# i=1118:1194 j=999
# V-component
tmp = df.thetao[:,:,998,1117:1194].compute()
tmp2 = tmp.groupby('time.year').mean('time')

# tmp=df.thetao[:,:,998,1117:1194]
# plt.pcolormesh(tmp.lon,-tmp.lev,tmp.data[0,:,:],cmap='needJet2');plt.colorbar();plt.ylim(-4850,0);



