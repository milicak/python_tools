import numpy as np 

gr = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')

root_folder = '/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_nemo4_2_01/rebuilt/'

ls1 = sorted(glob.glob(root_folder+'*grid_T.nc'))


fig, ax1 = plt.subplots(figsize=(12,8)) 
ind = 0
for ind3, fname in enumerate(ls1[120:155]):
    df = xr.open_dataset(fname)
    for ind2 in range(0,6):
        plt.pcolormesh(gr.nav_lon,gr.nav_lat,df.tos[ind2,:,:].where(df.tos[-1,:,:]!=0)
                       ,vmin=4,vmax=28);plt.colorbar();
        plt.title(str(ind))
        print(ind)
        ind += 1
        fout = 'gifs/BS_sst_' + str(ind).zfill(3) + '.png'          
        plt.savefig(fout, bbox_inches='tight',format='png',dpi=300) 
        plt.clf()
    
