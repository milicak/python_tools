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

def compute_dissipation_energy(uvel,vvel,wvel,u_rho,v_rho,w_rho,delta_x,delta_y,delta_z,visc):
# compute epsilon = tke dissipation
    # x-component
    dudx=(uvel[:,:,1:]-uvel[:,:,:-1])/delta_x;
    eps11=visc*(dudx**2)#*delta_x*delta_y*delta_z;

    dudy=np.zeros([nz,ny,nz])
    dudy[:,1:,:]=(uvel[:,1:,1:]-uvel[:,:-1,1:])/(0.5*(delta_y[:,1:,:]+delta_y[:,:-1,:]));
    dudy=np.concatenate((dudy, np.zeros([nz,ny,1])), axis=2)
    dudy[:,:,-1]=dudy[:,:,0]
    dudy=np.concatenate((dudy, np.zeros([nz,1,nx+1])), axis=1)
    dudy[:,-1,:]=dudy[:,0,:]
    dnm=visc*dudy**2;
    eps12=0.25*(dnm[:,:-1,:-1]+dnm[:,1:,1:]+dnm[:,1:,1:]+dnm[:,:-1,:-1])#*delta_x*delta_y*delta_z;

    dudz_w=np.zeros([nz,ny,nz])
    dudz_w[1:,:,:]=(u_rho[1:,:,:]-u_rho[:-1,:,:])/(0.5*(delta_z[1:,:,:]+delta_z[:-1,:,:]));
    dudz_w=np.concatenate((dudz_w, np.zeros([1,ny,nx])), axis=0)
    dudz=0.5*(dudz_w[:-1,:,:]+dudz_w[1:,:,:]);
    eps13=visc*(dudz**2)#*delta_x*delta_y*delta_z;

    # y-component

    dvdx=np.zeros([nz,ny,nz])
    dvdx[:,:,1:]=(vvel[:,1:,1:]-vvel[:,1:,:-1])/(0.5*(delta_x[:,:,1:]+delta_x[:,:,:-1]));
    dvdx=np.concatenate((dvdx, np.zeros([nz,1,nx])), axis=1)
    dvdx[:,-1,:]=dvdx[:,0,:]
    dvdx=np.concatenate((dvdx, np.zeros([nz,ny+1,1])), axis=2)
    dvdx[:,:,-1]=dvdx[:,:,0]
    dnm=visc*dvdx**2;
    eps21=0.25*(dnm[:,:-1,:-1]+dnm[:,1:,1:]+dnm[:,1:,1:]+dnm[:,:-1,:-1])#*delta_x*delta_y*delta_z;

    dvdy=(vvel[:,1:,:]-vvel[:,:-1,:])/delta_y;
    eps22=visc*(dvdy**2)#*delta_x*delta_y*delta_z;

    dvdz_w=np.zeros([nz,ny,nz])
    dvdz_w[1:,:,:]=(v_rho[1:,:,:]-v_rho[:-1,:,:])/(0.5*(delta_z[1:,:,:]+delta_z[:-1,:,:]));
    dvdz_w=np.concatenate((dvdz_w, np.zeros([1,ny,nx])), axis=0)
    dvdz=0.5*(dvdz_w[:-1,:,:]+dvdz_w[1:,:,:]);
    eps23=visc*(dvdz**2)#*delta_x*delta_y*delta_z;

    # z-component
    dwdx_u=np.zeros([nz,ny,nz])
    dwdx_u[:,:,1:]=(w_rho[:,:,1:]-w_rho[:,:,:-1])/(0.5*(delta_x[:,:,1:]+delta_x[:,:,:-1]));
    dwdx_u=np.concatenate((dwdx_u, np.zeros([nz,ny,1])), axis=2)
    dwdx=0.5*(dwdx_u[:,:,:-1]+dwdx_u[:,:,1:]);
    eps31=visc*(dwdx**2)#*delta_x*delta_y*delta_z;

    dwdy_v=np.zeros([nz,ny,nz])
    dwdy_v[:,1:,:]=(w_rho[:,1:,:]-w_rho[:,:-1,:])/(0.5*(delta_y[:,1:,:]+delta_y[:,:-1,:]));
    dwdy_v=np.concatenate((dwdy_v, np.zeros([nz,1,nx])), axis=1)
    dwdy=0.5*(dwdy_v[:,:-1,:]+dwdy_v[:,1:,:]);
    eps32=visc*(dwdy**2)#*delta_x*delta_y*delta_z;

    dwdz=(wvel[1:,:,:]-wvel[:-1,:,:])/delta_z;
    eps33=visc*(dwdz**2)#*delta_x*delta_y*delta_z;

    epsilon=eps11 + eps12 + eps13 + eps21 + eps22 + eps23 + eps31 +eps32 +eps33;
    #epsilon_rho=rho*(eps11 + eps12 + eps13 + eps21 + eps22 + eps23 + eps31 +eps32 +eps33);
    return epsilon

nx=512
ny=512
nz=512
cmap_needjet2=shfn()
visc=2e-5;
kappa=1.4e-5;
Pr=visc/kappa;
alpha_T=2; #alpha
beta_S=8; #beta
g=9.81;
rho0=1027;


root_folder='/export/grunchfs/unibjerknes/milicak/bckup/mitgcm/ice_leads/'

projects=['Exp01.3','Exp01.4','Exp01.5','Exp01.6','Exp01.7','Exp01.8','Exp01.9','Exp01.10','Exp01.11']

projectslbs=['Exp01_3','Exp01_4','Exp01_5','Exp01_6','Exp01_7','Exp01_8','Exp01_9','Exp01_10','Exp01_11']

itr=900*14
variable_name=['S']; #T for temp; S for salt

# IMPORTANT
# in MITgcm when you read it, the dimensions are Nz, Ny Nz

for i in range(0,9):
    print i,projects[i]
    foldername=root_folder+projects[i]+'/'
    if i==0:
        depth=rdmds(foldername+'Depth');
        xc=rdmds(foldername+'XC');
        yc=rdmds(foldername+'YC');
        dxc=rdmds(foldername+'DXC');
        dyc=rdmds(foldername+'DYC');
        drc=rdmds(foldername+'DRC');
        drf=rdmds(foldername+'DRF');
        Z=np.cumsum(drc);
        x=np.squeeze(xc[0,:])
        y=np.squeeze(yc[:,0])
        section=255
        delta_z=np.tile(np.array(drc),(ny,nx))
        delta_z_ref=np.copy(delta_z)
        delta_x=np.tile(dxc,(nz,1,1))
        delta_y=np.tile(dyc,(nz,1,1))

    salt = rdmds(foldername+'S',itr);
    uvel = rdmds(foldername+'U',itr);
    vvel = rdmds(foldername+'V',itr);
    wvel = rdmds(foldername+'W',itr);
    eta = rdmds(foldername+'Eta',itr);
    # This is due to zstar used in MITGCM; look at table 7.1 in mom4p1 manual
    corr = (1.0+eta/depth);
    #corr=repmat(corr,[1 1 Nz]);
    delta_z = delta_z_ref*corr;
    # for periodic and/or closed boundries in MITgcm
    uvel = np.concatenate((uvel, np.zeros([nz,ny,1])), axis=2)
    uvel[:,:,-1] = uvel[:,:,0]
    vvel = np.concatenate((vvel, np.zeros([nz,1,nx])), axis=1)
    vvel[:,-1,:] = vvel[:,0,:]
    wvel = np.concatenate((wvel, np.zeros([1,ny,nx])), axis=0)
    #wvel=np.append(wvel, np.zeros([1,ny,nx]), axis=0)
    u_rho = 0.5*(uvel[:,:,:-1]+uvel[:,:,1:]);
    v_rho = 0.5*(vvel[:,:-1,:]+vvel[:,1:,:]);
    w_rho = 0.5*(wvel[1:,:,:]+wvel[:-1,:,:]);

    # compute rho using linear EOS
    rho = (rho0-alpha_T*(-2)+beta_S*salt);

    epsilon = compute_dissipation_energy(uvel,vvel,wvel,u_rho,v_rho,w_rho,delta_x,delta_y,delta_z,visc)

    # xz section
    fig = plt.figure()
    #im1 = pcolor(x,-Z,np.squeeze(variable[:,section,:]),cmap=cmap_needjet2,vmin=32,vmax=32.02)
    im1 = plt.pcolormesh(x,-Z,np.log10(np.squeeze(epsilon[:,section,:])),linewidth=0,rasterized=True,shading='flat',cmap="jet",vmin=-10,vmax=-6)
    #im1.set_edgecolor('face')
    plt.ylim((-128,0))
    plt.xlim((0,128))
    cb = plt.colorbar(im1,pad=0.02) # pad is the distance between colorbar and figure
    cb.set_label(r'log$_{10} \varepsilon$' ' [W/kg]')
#    cb.set_label('[' r'$^\circ$' 'C]')
    plt.ylabel('depth [m]')
    plt.xlabel('x [m]')
    #plt.show()
    plt.savefig('paperfigs/verticalxz_epsilon_section_'+projectslbs[i]+'_'+str(itr)+'.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)
    # yz section
    fig = plt.figure()
    #im1 = pcolor(x,-Z,np.squeeze(variable[:,section,:]),cmap=cmap_needjet2,vmin=32,vmax=32.02)
    im1 = plt.pcolormesh(y,-Z,np.log10(np.squeeze(epsilon[:,:,section])),linewidth=0,rasterized=True,shading='flat',cmap="jet",vmin=-10,vmax=-6)
    #im1.set_edgecolor('face')
    plt.ylim((-128,0))
    plt.xlim((0,128))
    cb = plt.colorbar(im1,pad=0.02) # pad is the distance between colorbar and figure
    cb.set_label(r'log$_{10} \varepsilon$' ' [W/kg]')
#    cb.set_label('[' r'$^\circ$' 'C]')
    plt.ylabel('depth [m]')
    plt.xlabel('y [m]')
    #plt.show()
    plt.savefig('paperfigs/verticalyz_epsilon_section_'+projectslbs[i]+'_'+str(itr)+'.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)


