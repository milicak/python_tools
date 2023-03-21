import numpy as np 

gr = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')

# root_folder = '/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/rebuilt/'
root_folder = '/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_bdy7_geom3_rlx/rebuilt/'
root_folder = '/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_bdy7_geom3_rlxbiharm/rebuilt/'

ls1 = sorted(glob.glob(root_folder+'*grid_T.nc'))

ds = xr.open_dataset('BlackSea_mask.nc')
df = xr.open_dataset(ls1[0])
mask = xr.where(df.thetao[0,:,:,:]!=0,1,0)
mask = mask*ds.mask
area = gr.e1t*gr.e2t

fig, ax1 = plt.subplots(figsize=(8,8)) 
for ind, fname in enumerate(ls1):
    df = xr.open_dataset(fname)
    print(ind)
    plt.pcolormesh(gr.nav_lon,gr.nav_lat,df.so[-1,55,:,:].where(mask[55,:,:]!=0),vmin=17,vmax=21);plt.colorbar();
    plt.title(str(ind))
    fout = 'gifs/BS_salt_' + str(ind).zfill(3) + '.png'          
    plt.savefig(fout, bbox_inches='tight',format='png',dpi=300) 
    plt.clf()

