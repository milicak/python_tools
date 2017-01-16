''' Create an N-bin discrete colormap from the specified input map '''
import matplotlib.pyplot as plt
import numpy as np


def discrete_cmap(Nnum, base_cmap=None):
    """Create an Nnum-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, Nnum))
    cmap_name = base.name + str(Nnum)
    return base.from_list(cmap_name, color_list, Nnum)
