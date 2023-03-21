import os

intermediate2 = xr.open_dataset('bs-test-int-bdy7_argo_2018.nc')
intermediate = intermediate2.where((intermediate2['latitude'] > 42.) &
                                  (intermediate2['latitude'] < 44.) &
                                  (intermediate2['longitude'] > 35.) &
                                  (intermediate2['longitude'] < 38.5)).dropna('obs')
intermediate
>
        (obs: 31242, model: 1)
        (model) object 'bs-test-int-bdy7'
        (obs) float32 36.9 36.9 36.9 36.9 ... 37.99 37.99 37.99
        (obs) float64 7.738 13.0 17.46 ... 967.0 972.0 976.0
        (obs) datetime64[ns] 2018-01-06T00:05:00 ... 2018-12-2...
e       (obs) int64 54599821 54599821 ... 62827631 62827631
        (obs) float32 42.98 42.98 42.98 ... 43.62 43.62 43.62
out coordinates: obs

        (obs) float64 11.86 11.86 11.86 ... 8.953 8.954 8.954
        (obs) float64 18.39 18.39 18.39 ... 22.27 22.27 22.27
rature  (model, obs) float64 10.52 10.52 10.52 ... 8.936 8.936
ity     (model, obs) float64 19.01 19.01 19.01 ... 22.26 22.26
for dc_reference, profile in intermediate.groupby('dc_reference'):
    if profile['time.season'][0]=='DJF':
        plt.plot(profile.salinity,-profile.depth,color='yellow')


