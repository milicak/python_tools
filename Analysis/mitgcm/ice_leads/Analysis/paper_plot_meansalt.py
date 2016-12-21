import numpy as np
#%matplotlib inline
#np.shape !!!!!

import matplotlib.pyplot as plt
import scipy.io
import numpy.ma as ma
#import my_nanfilter
#from my_nanfilter import my_nanfilterbox
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
# MITGCM packages
sys.path.append('/home/mil021/models/MITgcm/utils/MITgcmutils/')
from MITgcmutils import rdmds
from needJet2 import shfn
from disc_cb import discrete_cmap




#def runningMeanFast(x, N):
#    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def plt_arrays(x, y, title="", label="", color="red", linestyle="solid", linewidth=2):    
    axes = fig.add_subplot(111)
    axes.plot(x,y, color=color, linestyle=linestyle, linewidth=linewidth,label=label)
    axes.set_title(title)
    #axes.set_label(label)
    #plt.legend(label)
    axes.grid()
    legend = axes.legend(loc='upper right')


fyear = '1'; # first year
lyear = '250'; # last year
dn = 5; # 5 year filter
time = np.linspace(1, 250, 250)


cmap_needjet2=shfn()

root_folder='/mnt/fimmhome/Analysis/mitgcm/ice_leads/Analysis/matfiles/'

#projects=['Exp01.3','Exp01.4','Exp01.5','Exp01.6','Exp01.7','Exp01.8','Exp01.9','Exp01.10','Exp01.11']
#projectslbs=['Exp01_3','Exp01_4','Exp01_5','Exp01_6','Exp01_7','Exp01_8','Exp01_9','Exp01_10','Exp01_11']

projects=['Exp01.3','Exp01.4','Exp01.5','Exp01.9','Exp01.10','Exp01.11']
projectslbs=['Exp01_3','Exp01_4','Exp01_5','Exp01_9','Exp01_10','Exp01_11']
plotcolors=['cyan','blue','magenta','black','green','red','brown']

Z=np.linspace(0.0,128.0,512)
            
fig = plt.figure()
for i in xrange(0,6):
    filename=root_folder+projects[i]+'_mean_salt.mat'
    mat = scipy.io.loadmat(filename)
    meansalt=np.transpose(np.array(mat['meansalt']))    
    plt_arrays(meansalt[:,-1],-Z, title="",label=projectslbs[i],color=plotcolors[i])
plt.xlabel('Salinity [psu]')    
plt.ylabel('depth [m]')
plt.xticks([32,32.005,32.01,32.015,32.02,32.025],('32','32.005','32.01','32.015','32.02','32.025'))
plt_arrays(meansalt[:,0],-Z, title="",color="grey",linestyle="dashed")
plt.legend(loc='lower left') 
plt.show()
plt.savefig('paperfigs/meansalt.eps', format='eps', dpi=300)    
plt.close(fig)


fig = plt.figure()
for i in xrange(0,6):
    filename=root_folder+projects[i]+'_energetics_dimensional.mat'
    mat = scipy.io.loadmat(filename)
    rpe=np.transpose(np.array(mat['rpe']))    
    time_days=np.transpose(np.array(mat['time_days']))    
    plt_arrays(time_days,(rpe-rpe[0])/rpe[0], title="",label=projectslbs[i],color=plotcolors[i])
plt.xlabel('time [days]')    
plt.ylabel('(RPE-RPE(0))/RPE(0)')
#plt.xticks([32,32.005,32.01,32.015,32.02,32.025],('32','32.005','32.01','32.015','32.02','32.025'))
#plt_arrays(meansalt[:,0],-Z, title="",color="grey",linestyle="dashed")
plt.legend(loc='upper left') 
plt.xlim(0,1.6)
plt.show()
plt.savefig('paperfigs/rpe_time.eps', format='eps', dpi=300)    
plt.close(fig)