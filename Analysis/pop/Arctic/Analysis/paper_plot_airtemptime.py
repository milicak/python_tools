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

#root_folder='/mnt/fimm/Analysis/pop/Arctic/Analysis/matfiles/'
root_folder='/export/grunchfs/unibjerknes/milicak/bckup/Analysis/pop/Arctic/Analysis/matfiles/'

projects=['B1850CN_f19_tn11_kdsens','B1850CN_f19_tn11_kdsens01','B1850CN_f19_tn11_kdsens02',
          'B1850CN_f19_tn11_kdsens03','B1850CN_f19_tn11_kdsens05','B1850CN_f19_tn11_kdsens04',
          'B1850CN_f19_tn11_kdsens06','B1850CN_f19_tn11_kdsens07']

legendnames=['Cold-ctrl','Exp1','Exp2','Exp3','Warm-ctrl','Exp4','Exp5','Exp6']
plotcolors=['cyan','blue','magenta','brown','green','red','black','orange']


fig = plt.figure()
for i in xrange(0,4):
    filename=root_folder+projects[i]+'_airtemp_mean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    AIRTEMPTrefwinter=np.array(mat['AIRTEMPTrefwinter'])
    dnm=np.squeeze(AIRTEMPTrefwinter[:,:,48:])
    dnm=np.squeeze(np.mean(dnm,2))
    dnm=np.squeeze(np.mean(dnm,1))
    AIRTEMPTrefwinter_fil=my_nanfilterbox(dnm,dn)
    #plt_arrays(time,runningMeanFast(amoc26[:,0],10), title="",label=legendnames[i],color=plotcolors[i])
    plt_arrays(time,AIRTEMPTrefwinter_fil, title="",label=legendnames[i],color=plotcolors[i])
plt.xlabel('Time [years]')    
plt.ylabel('Northern Hemisphere air temperature [' r'$^\circ$' 'C]') 
#plt.ylim(17,24) 
plt.show()
plt.savefig('paperfigs/cold_exps_airtempwinternorth_time.eps', format='eps', dpi=1000)    
plt.close(fig)

fig = plt.figure()
for i in xrange(4,8):
    filename=root_folder+projects[i]+'_airtemp_mean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    AIRTEMPTrefwinter=np.array(mat['AIRTEMPTrefwinter'])
    dnm=np.squeeze(AIRTEMPTrefwinter[:,:,48:])
    dnm=np.squeeze(np.mean(dnm,2))
    dnm=np.squeeze(np.mean(dnm,1))
    AIRTEMPTrefwinter_fil=my_nanfilterbox(dnm,dn)
    #plt_arrays(time,runningMeanFast(amoc26[:,0],10), title="",label=legendnames[i],color=plotcolors[i])
    plt_arrays(time,AIRTEMPTrefwinter_fil, title="",label=legendnames[i],color=plotcolors[i])
plt.xlabel('Time [years]')    
plt.ylabel('Northern Hemisphere air temperature [' r'$^\circ$' 'C]')  
#plt.ylim(17,24) 
plt.show()
plt.savefig('paperfigs/warm_exps_airtempwinternorth_time.eps', format='eps', dpi=1000)    
plt.close(fig)

ii=[0,4]
fig = plt.figure()
for i in ii:
    print i
    filename=root_folder+projects[i]+'_airtemp_mean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    AIRTEMPTrefwinter=np.array(mat['AIRTEMPTrefwinter'])
    dnm=np.squeeze(AIRTEMPTrefwinter[:,:,48:])
    dnm=np.squeeze(np.mean(dnm,2))
    dnm=np.squeeze(np.mean(dnm,1))
    AIRTEMPTrefwinter_fil=my_nanfilterbox(dnm,dn)
    #plt_arrays(time,runningMeanFast(amoc26[:,0],10), title="",label=legendnames[i],color=plotcolors[i])
    plt_arrays(time,AIRTEMPTrefwinter_fil, title="",label=legendnames[i],color=plotcolors[i])
plt.xlabel('Time [years]')    
plt.ylabel('Northern Hemisphere air temperature [' r'$^\circ$' 'C]')  
#plt.ylim(17,24) 
plt.show()
plt.savefig('paperfigs/cold_warm_ctrls_airtempwinternorth_time.eps', format='eps', dpi=1000)    
plt.close(fig)
