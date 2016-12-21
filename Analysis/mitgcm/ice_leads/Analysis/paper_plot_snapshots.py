import numpy as np
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.io
import numpy.ma as ma
from disc_cb import discrete_cmap
#import my_nanfilter
#from my_nanfilter import my_nanfilterbox
import nccf
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

# MITGCM packages
sys.path.append('/fimm/home/bjerknes/milicak/models/MITgcm/utils/MITgcmutils/')
from MITgcmutils import rdmds
from needJet2 import shfn
from disc_cb import discrete_cmap

nx=512
ny=512
nz=512
cmap_needjet2=shfn()

root_folder='/export/grunchfs/unibjerknes/milicak/bckup/mitgcm/ice_leads/'

projects=['Exp01.3','Exp01.4','Exp01.5','Exp01.6','Exp01.7','Exp01.8','Exp01.9','Exp01.10','Exp01.11']

projectslbs=['Exp01_3','Exp01_4','Exp01_5','Exp01_6','Exp01_7','Exp01_8','Exp01_9','Exp01_10','Exp01_11']

itr=900*14
variable_name=['S']; #T for temp; S for salt

# compute amoc    
for i in range(0,9):
    print i,projects[i]
    foldername=root_folder+projects[i]+'/'
    print foldername
    if i==0:
        depth=rdmds(foldername+'Depth');
        xc=rdmds(foldername+'XC');
        yc=rdmds(foldername+'YC');
        drc=rdmds(foldername+'DRC');
        Z=np.cumsum(drc);
        x=np.squeeze(xc[0,:])
        y=np.squeeze(yc[:,0])
        section=255
    
    
    variable=rdmds(foldername+'S',itr);
    # xz section
    fig = plt.figure()
    #im1 = pcolor(x,-Z,np.squeeze(variable[:,section,:]),cmap=cmap_needjet2,vmin=32,vmax=32.02)
    im1 = plt.pcolormesh(x,-Z,np.squeeze(variable[:,section,:]),linewidth=0,rasterized=True,shading='flat',cmap=cmap_needjet2,vmin=32,vmax=32.02)
    #im1.set_edgecolor('face')
    plt.ylim((-128,0))           
    plt.xlim((0,128))    
    cb = plt.colorbar(im1,pad=0.02) # pad is the distance between colorbar and figure
    cb.set_label('[psu]')
#    cb.set_label('[' r'$^\circ$' 'C]')
    plt.ylabel('depth [m]')    
    plt.xlabel('x [m]')
    #plt.show() 
    plt.savefig('paperfigs/verticalxz_section_'+projectslbs[i]+'_'+str(itr)+'.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)
    # yz section
    fig = plt.figure()
    #im1 = pcolor(x,-Z,np.squeeze(variable[:,section,:]),cmap=cmap_needjet2,vmin=32,vmax=32.02)
    im1 = plt.pcolormesh(y,-Z,np.squeeze(variable[:,:,section]),linewidth=0,rasterized=True,shading='flat',cmap=cmap_needjet2,vmin=32,vmax=32.02)
    #im1.set_edgecolor('face')
    plt.ylim((-128,0))           
    plt.xlim((0,128))    
    cb = plt.colorbar(im1,pad=0.02) # pad is the distance between colorbar and figure
    cb.set_label('[psu]')
#    cb.set_label('[' r'$^\circ$' 'C]')
    plt.ylabel('depth [m]')    
    plt.xlabel('y [m]')
    #plt.show() 
    plt.savefig('paperfigs/verticalyz_section_'+projectslbs[i]+'_'+str(itr)+'.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)