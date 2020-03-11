"""
remove average from sea level 
at boundary input of SHYFEM

this is specific for uTSS
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
plt.ion()

# make plots
make_plots = False
# make_plots = True

# mnthnames = ['20160101', '20160201', '20160301', '20160401',
             # '20160501', '20160601', '20160701', '20160801',
             # '20160901', '20161001', '20161101', '20161201']
# mnthnames = ['20162017years_analysis']
mnthnames = ['201620172018years_analysis']
# mnthnames = ['dnm']
# mnthnames = ['2016year']


for foldername in mnthnames:
    pathname = "/work/mi19918/Projects/uTSS/preproc/runs/"+foldername
    print 'directory '+foldername
    os.chdir(pathname)
    
    header 		= 3	# number of header lines in boundary files
    boundaries 	= 3	# number of boundary files to handle
    
    # date of first record
    date0		= datetime(2016,1,1,0,0,0)
    delta		= timedelta(seconds=86400)
    
    # filenames and number of boundary nodes of
    # each boundary file
    fbound1, nbound1 	= 'boundn_L1_TSS.dat', 26
    fbound2, nbound2 	= 'boundn_L2_TSS.dat', 55 
    fbound3, nbound3 	= 'boundn_L3_TSS.dat', 92
    
    # name of new files
    fbound1a = fbound1.replace('_TSS.dat','anew_TSS.dat')
    fbound2a = fbound2.replace('_TSS.dat','anew_TSS.dat')
    fbound3a = fbound3.replace('_TSS.dat','anew_TSS.dat')
    
    newfiles = [fbound1a,fbound2a,fbound3a]
    
    nbounds = [nbound1,nbound2,nbound3]
    
    with open(fbound1,'r') as f1: lines1 = f1.readlines()
    with open(fbound2,'r') as f2: lines2 = f2.readlines()
    with open(fbound3,'r') as f3: lines3 = f3.readlines()
    
    liness = (lines1,lines2,lines3)
    
    ## lines for each time step (header + values ), is different for each boundary file
    lin1 = header + nbound1
    lin2 = header + nbound2
    lin3 = header + nbound3
    
    ## times steps included in file
    nt1 = len(lines1) / lin1
    nt2 = len(lines2) / lin2
    nt3 = len(lines3) / lin3
    
    print '%s time steps' % nt1
    print '%s time steps' % nt2
    print '%s time steps' % nt3
    
    nts = [nt1,nt2,nt3]
    
    # check that all nts are equal
    check = nts.count(nts[0]) == len(nts)
    
    print check
    
    if check == False:
    	print 'time steps of boundary files are not the same'
    	print 'check time of boundary files'
    	sys.exit()
    
    # create list of dates 
    dates = [date0+delta*i for i in range(nt1)]
    
    # store boundaries in numpy arrays
    boundn1 	= np.zeros((nt1,nbound1))
    boundn2 	= np.zeros((nt1,nbound2))
    boundn3 	= np.zeros((nt1,nbound3))
    
    for bb in range(boundaries):
    	print 'boundary L%s: ', bb+1
    	for t in range(nt1):
    	        i1 = t*(header+nbounds[bb])+header
           	i2 = i1 + nbounds[bb]
            	vals = liness[bb][i1:i2]
    		if bb == 0: boundn1[t,:] = np.asarray([float(i) for i in vals])
    		elif bb == 1: boundn2[t,:] = np.asarray([float(i) for i in vals])
    		elif bb == 2: boundn3[t,:] = np.asarray([float(i) for i in vals])
    
    # calculate average over period of L1 and L2 (weighted average)
    w1,w2 		= float(nbound1)/float(nbound1+nbound2),float(nbound2)/float(nbound1+nbound2)
    l1l2_ave 	= w1 * np.average(boundn1) + w2 * np.average(boundn2)
    l1l2_min 	= w1 * np.min(boundn1) + w2 * np.min(boundn2)
    
    print 'l1l2_ave : %s ' % l1l2_ave, l1l2_min
    # print 'mehmet change with min '
    # l1l2_ave = l1l2_min
    # print 'l1l2_ave : %s ' % l1l2_ave
    # remove average Black Sea mean ssh from boundary values
    # 2.06 was the 2016 annual mean of Basin BS ssh 
    # boundn3    = boundn3 - 2.06 
    # 2016 and 2017 working one is 2.32
    boundn3[:365+366,:]    = boundn3[:365+366,:] - 2.32 
    diff = boundn3[365+366-1]-boundn3[365+366]
    boundn3[365+366:]+=diff
    
    # remove average from boundaries
    boundn1a	= boundn1 - l1l2_ave
    boundn2a	= boundn2 - l1l2_ave
    boundn3a	= boundn3 - l1l2_ave
    # boundn1a	= boundn1 - l1l2_min
    # boundn2a	= boundn2 - l1l2_min
    # boundn3a	= boundn3 - l1l2_min
    #boundn1a	= boundn1 + 0.5 
    #boundn2a	= boundn2 + 0.5
    adddlin = (np.arange(550,boundn3.shape[0],1)-550)*0.15/(boundn3.shape[0]-550)
    aa = np.tile(adddlin,(92,1))
    boundn3a[550::,:] = boundn3a[550::,:]+0.3+np.transpose(aa)
        
    
    # coordinates for boundaries in plot
    x1 = np.arange(nbound1)
    x2 = np.arange(nbound1,nbound1+nbound2)
    x3 = np.arange(nbound1+nbound2,nbound1+nbound2+nbound3)
    
    
    if make_plots == True:
    	print 'making plots of boundaries'
    	for tt in range(nt1):
    		fig,ax 	= plt.subplots(figsize=(10,7))
    		l1, 	= ax.plot(x1,boundn1[tt,:])
    		l2, 	= ax.plot(x2,boundn2[tt,:])
    		l3, 	= ax.plot(x3,boundn3[tt,:])
    		l1a, 	= ax.plot(x1,boundn1a[tt,:])
    		l2a, 	= ax.plot(x2,boundn2a[tt,:])
    		l3a, 	= ax.plot(x3,boundn3a[tt,:])
    
    		#ax.legend([l1,l2,l3],['L1','L2','L3'],loc='upper left')
    		ax.legend([l1,l2,l3,l1a,l2a,l3a],['L1','L2','L3','L1a','L2a','L3a'],loc='upper left')
    		plt.title('boundaries at %s' % dates[tt].strftime('%Y/%m/%d - %H:%M:%S'))
    		plt.xlabel('grid points')
    		plt.ylabel('[m]')
    		plt.ylim([-1,1])
    		plt.savefig('bound_plots/boundaries_t_%s.png' % tt )
    		plt.close()
    
    #for t in range(nt):
    #        i1 = t*(header+nbound)+header
    #        i2 = i1 + nbound
    #        vals = lines[i1:i2]
    #        ave = np.mean(np.asarray([float(i) for i in vals]))
    #	print lines[i1-2],ave
    #        vals = [str(float(i)-ave)+'\n' for i in vals]
    #        lines[i1:i2] = vals
    
    for bb in range(boundaries):
    	print 'shifting boundary L%s ...', bb+1
    	for t in range(nt1):
    	        i1 = t*(header+nbounds[bb])+header
           	i2 = i1 + nbounds[bb]
    		if bb == 0: vals = [str(i)+'\n' for i in boundn1a[t,:]]
    		elif bb == 1: vals = [str(i)+'\n' for i in boundn2a[t,:]]
    		elif bb == 2: vals = [str(i)+'\n' for i in boundn3a[t,:]]
    		liness[bb][i1:i2] = vals		
    
    print 'newfiles ', newfiles
    
    for bb in range(boundaries):
    	print 'writing new boundary file for L%s ', bb+1
    	print 'new filename: ', newfiles[bb]
    	g = open(newfiles[bb],'w')
    	for i in range(len(liness[bb])):
    		g.write(liness[bb][i])
    	g.close() 
    
    




