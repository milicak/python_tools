import numpy as np
import matplotlib as mpl
def shfn():
      ff = np.loadtxt('/fimm/home/bjerknes/milicak/matlab/tools/jet4_python')
      cmap_needjet2= mpl.colors.LinearSegmentedColormap.from_list("my_colormap",ff, N=54, gamma=1.0)
      return cmap_needjet2

