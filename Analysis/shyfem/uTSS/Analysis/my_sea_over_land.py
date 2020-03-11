import numpy as np

def seaoverland_2d(data, iterations=1, copy=False):

    if copy:
        data = np.ma.copy(data)

    if not isinstance(data, np.ma.masked_array) or not data.mask.any():
        return data

    for _ in range(iterations):
        shifted = []
        ni, nj = data.shape
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i != 0 or j != 0:
                    # Shift the grid by i horizontally and j vertically and
                    # append it to the array. Shifted grids are 2 units smaller
                    # in both dimensions to accomodate the shift.
                    shifted.append(data[1 + i:ni - 1 + i, 1 + j:nj - 1 + j])

        # Calculate the mean value of the shifted grids to obtain the
        # approximated values. Only non-masked entries are taken into account.
        approx = np.ma.mean(shifted, axis=0)

        # Create a view without the outer points (so it is the same size as the
        # shifted grids), then copy the approximated values for the cells that
        # are masked.
        view = data[1:-1, 1:-1]
        np.copyto(view, approx, where=(view.mask & ~approx.mask))

        # Combine the two masks, unmasking values that were masked in view but
        # have been successfully approximated.
        view.mask &= approx.mask

    #brutale ma necessario
    data[0,:] = data[1,:]
    data[-1,:] = data[-2,:]
    data[:,0] = data[:,1]
    data[:,-1] = data[:,-2]


    return data


def seaoverland_1d(data, iterations=1, copy=False):

    if copy:
        data = np.ma.copy(data)

    if not isinstance(data, np.ma.masked_array) or not data.mask.any():
        return data

    for _ in range(iterations):
        shifted = []
        nj = len(data)
        for j in range(-1, 2):
            	if j != 0:
                    	# Shift the grid by horizontally and
                    	# append it to the array. Shifted grids are 2 units smaller
                    	# in both dimensions to accomodate the shift.
                    	shifted.append(data[ 1 + j:nj - 1 + j])

        # Calculate the mean value of the shifted grids to obtain the
        # approximated values. Only non-masked entries are taken into account.
        approx = np.ma.mean(shifted, axis=0)

        # Create a view without the outer points (so it is the same size as the
        # shifted grids), then copy the approximated values for the cells that
        # are masked.
        view = data[ 1:-1]
        np.copyto(view, approx, where=(view.mask & ~approx.mask))

        # Combine the two masks, unmasking values that were masked in view but
        # have been successfully approximated.
        view.mask &= approx.mask

    # brutal but necessary		
    data[0] = data[1]
    data[-1] = data[-2]

    return data

