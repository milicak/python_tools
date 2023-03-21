import xarray as xr
import numpy as np
from pqtool.utils import metrics
import matplotlib.pyplot as plt

root_folder = '/work/opa/mi19918/intermediate/'
# expid = 'BS-NRT_MI_bdy7_geom3_rlxbiharm02'
# expid = 'bs-test-int-bdy7'
expid = 'BS-BIHARM_003'
bins = np.array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160,
                 180, 200, 300, 400, 500, 1000])
yr = 2014
tind = 1
for year in range(yr,yr+5):
    intermediate = xr.open_dataset(root_folder + expid + '_argo_' +
                                     str(year) + '.nc')

    result = intermediate.groupby_bins('depth', bins=bins).apply(metrics)
    # result.to_dataframe()
    for variable in ['temperature']:
    # for variable in ['salinity']:
    # for variable in ['temperature', 'salinity']:
        for metric in ['Bias']:
        # for metric in ['Bias', 'RMSE']:
            if variable == 'temperature' and metric == 'Bias':
                figno = 1
                plt.figure(figno)
                if tind == 1:
                    ax1 = plt.axes()
                    tind += 1
            if variable == 'temperature' and metric == 'RMSE':
                figno = 2
                fig = plt.figure(figno)
                ax2 = plt.axes()
            if variable == 'salinity' and metric == 'Bias':
                figno = 3
                fig = plt.figure(figno)
                ax3 = plt.axes()
            if variable == 'salinity' and metric == 'RMSE':
                figno = 4
                fig = plt.figure(figno)
                ax4 = plt.axes()

            result['depth'] = xr.DataArray((bins[1:] + bins[:-1])/ 2., dims='depth_bins')
            ax1.plot(model, result['depth'], '.-', label=str(year))
            # for model in result['%s_%s' % (variable, metric.lower())].transpose('model', 'depth_bins'):
            #     print(year)
            #     ax1.plot(model, result['depth'], '.-', label=str(year))


plt.title(variable.capitalize())
plt.xlabel('%s [%s]' % (metric, intermediate[variable].attrs['units']))
plt.ylabel('Depth [m]')
plt.semilogy()
ax1.invert_yaxis()
ax1.set_yticks(bins)
ax1.set_yticklabels(bins)
plt.legend(loc='best')
plt.grid(True, c='silver', lw=1, ls=':')


# for figno in range(1,5):
#     fig = plt.figure(figno)
#     ax = plt.axes()
#     if figno<3:
#         variable = 'temperature'
#     else:
#         variable = 'salinity'
#
#     if figno == 1:
#         metric = 'Bias'
#     if figno == 2:
#         metric = 'RMSE'
#     if figno == 3:
#         metric = 'Bias'
#     if figno == 4:
#         metric = 'RMSE'
#
#     plt.title(variable.capitalize())
#     plt.xlabel('%s [%s]' % (metric, intermediate[variable].attrs['units']))
#     plt.ylabel('Depth [m]')
#     plt.semilogy()
#     ax.invert_yaxis()
#     ax.set_yticks(bins)
#     ax.set_yticklabels(bins)
#     plt.legend(loc='best')
#     plt.grid(True, c='silver', lw=1, ls=':')
