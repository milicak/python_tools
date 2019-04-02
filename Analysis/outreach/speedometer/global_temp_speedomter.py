import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
from numpy import loadtxt
plt.ioff()

def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1



    start_ind=np.int32(xdata.shape)-2
    end_ind=np.int32(xdata.shape)-1
    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color,linewidth=10.0),
        size=size,
    )



names=['','','1.7','1.5','1.3','1.1','0.9','0.7','0.5','0.3','0.1','','']
size=[10, 10, 10, 10, 10, 10,10, 10, 10, 10, 10, 55, 55]
ds = loadtxt("647_Global_Temperature_Data_File.txt")
ds2 = loadtxt('RCP85_annual_rel_1850_1900')
year = ds[:,0]
Temp = ds[:,1]

year = np.concatenate((ds[:,0],ds2[12:,0]))
Temp = np.concatenate((ds[:,1],ds2[12:,1]))

# angle is between -1 and 2
angles = 0.5*(2.0-Temp)*np.pi
angles = 0.5*(np.maximum(2.-Temp,0))*np.pi

colors=[(255./255, 0, 0),(244./255, 122./255, 66./255),(255./255, 138./255, 71./255),
       (255./255, 255./255, 71./255), (230./255, 255./255, 71./255),
       (166./255, 255./255, 71./255), (89./255, 255./255, 71./255),
        (71./255, 233./255, 255./255), (71./255, 218./255, 255./255),
        (71./255, 187./255, 255./255), (71./255, 132./255, 255./255), (159./255, 185./255, 242./255),'white']

ind = 1
for angle in angles:
    plt.figure(figsize=(10,10))
    plt.title('Global Annual Mean Surface Air Temperature Change, Year ='
               +np.str(year[ind-1]),size=18)
    t = np.linspace(0.0, .56, 200)
    x = np.cos(angle)*t
    y = np.sin(angle)*t
    line = plt.plot(x, y,linewidth=10.0,color='k')[0]
    #add_arrow(line,position=np.max(x)-0.001,size=150)
    add_arrow(line,position=np.min(x)+0.001,size=150)
    my_circle=plt.Circle( (0,0), 0.4, color='white')
    wedges,texts=plt.pie(np.transpose(size), labels=names, colors=colors,textprops=dict(color="k"),labeldistance=0.70);
    plt.text(0.5,0.1,'Are you insane!',size=18,rotation=8)
    plt.text(0.53,0.31,'Hot! Hot!',size=18,rotation=24)
    plt.text(-0.38,-0.0,'North of the Wall',size=18,rotation=-52)
    plt.setp(texts,size=18);
    plt.axis('equal');
    plt.tight_layout();
    p=plt.gcf()
    p.gca().add_artist(my_circle);
    if year[ind-1] <= 2017.0:
        plt.text(-0.25,1.03,'NASA GISTEMP',size=18)
    else:
        plt.text(-0.25,1.03,'NORESM RCP8.5',size=18)


    sdate = "%4.4d" % (ind)
    savename = 'gifs/global_temp_'+sdate+'.png'
    plt.savefig(savename,dpi=200,bbox_inches='tight',format='png')
    plt.close()
    print(ind)
    ind += 1



