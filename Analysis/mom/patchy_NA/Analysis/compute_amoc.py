import numpy as np
#%matplotlib inline
#np.shape !!!!!

import matplotlib.pyplot as plt
import scipy.io
import numpy.ma as ma
import my_nanfilter
from my_nanfilter import my_nanfilterbox
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

reload(my_nanfilter)

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

mask_file='/mnt/fimmhome/Analysis/mom/APE/SO/Analysis/grid_spec_v6_regMask.nc';
mask=nc_read(mask_file,'tmask');
# Atlantic mask==2

root_folder='/mnt/fimm/mom/'

projects=['om3_core3_ctrl','om3_core3_patchy_full_01','om3_core3_patchy_full_02']

legendnames=['Cold-ctrl','Exp1','Exp2','Exp3','Warm-ctrl','Exp4','Exp5']
plotcolors=['cyan','blue','magenta','brown','green','red','black']
hist_folder = ['history_63-124years'];
hist_folder0 = ['history_1-62years'];

# compute amoc    
for i in range(0,1):
    filename0=root_folder+projects[i]+'/'+hist_folder0[0]+'/'+'00010101.ocean_month.nc'
    filename=root_folder+projects[i]+'/'+hist_folder[0]+'/'+'00630101.ocean_month.nc'            
    tmp=nc_read(filename0,'ty_trans')
    amoc=np.copy(tmp)
    tmp=nc_read(filename,'ty_trans')    
    