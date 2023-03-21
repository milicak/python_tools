import xarray as xr
import numpy as np
from pqtool.utils import metrics
import matplotlib.pyplot as plt

root_folder = '/work/opa/mi19918/intermediate/'

yr = 2018
for year in range(yr,yr+1):
    # intermediate = xr.open_mfdataset(root_folder + '*' + np.str(year) + '.nc', combine='by_coords')
    intermediate = xr.open_mfdataset(root_folder + '*_argo_' + str(year) + '.nc', combine='by_coords')

for index, model in enumerate(intermediate['model']):
    print('%d: %s' % (index, model.data))

# intermediate = intermediate.isel(model=[3,4,5])    
# intermediate = intermediate.isel(model=[6,9, 10, 11, 12])    
# intermediate = intermediate.isel(model=[6,12,14,15,16,17])    
# intermediate = intermediate.isel(model=[17,18,19])
intermediate = intermediate.isel(model=[15,16,17])    
# Drop observations if they are not present for all models
intermediate = intermediate.dropna('obs')

bins = np.array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160,
                 180, 200, 300, 400, 500, 1000])
result = intermediate.groupby_bins('depth', bins=bins).apply(metrics)
result.to_dataframe()

for variable in ['temperature', 'salinity']:
    for metric in ['Bias', 'RMSE']:
        fig = plt.figure(dpi=120)
        ax = plt.axes()

        result['depth'] = xr.DataArray((bins[1:] + bins[:-1])/ 2., dims='depth_bins')
        for model in result['%s_%s' % (variable, metric.lower())].transpose('model', 'depth_bins'):
            plt.plot(model, result['depth'], '.-', label=model['model'].data)

        plt.semilogy()
        ax.invert_yaxis()
        ax.set_yticks(bins)
        ax.set_yticklabels(bins)

        plt.title(variable.capitalize())
        plt.xlabel('%s [%s]' % (metric, intermediate[variable].attrs['units']))
        plt.ylabel('Depth [m]')
        plt.legend(loc='best')
        plt.grid(True, c='silver', lw=1, ls=':')


