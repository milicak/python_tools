import numpy as np
def my_nanfilterbox(y, dn):
	# box filter
	filtermi = np.ones(dn)
	# bar filter
	#filtermi = np.barlett(dn)
	ny=np.max(np.size(y))-1
	dummy=np.NaN*np.ones(ny+dn)
	if np.float(dn)/2 == np.round(np.float(dn)/2):
	   dn0 = dn/2-1
  	else:
           dn0 = (dn-1)/2
      
        dummy[0+dn0:ny+dn0+1] = y
        yfilter=np.ones(ny+1)
        for n in range(0,ny+1):
	  fy = dummy[0+n:dn+n]*filtermi
	  II = (~np.isnan(fy))
 	  wt0 = 1/np.mean(filtermi[II])	
          yfilter[n] = wt0*np.mean(fy[II])
        
	return yfilter
