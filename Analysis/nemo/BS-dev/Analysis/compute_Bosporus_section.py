import numpy as np

y1 = np.arange(8,22)
x1 = 70*np.ones(y1.shape)
y2 = np.array([21, 22, 22, 23])
x2 = np.array([71, 71, 72, 72])
y3 = np.arange(23,28)
x3 = 73*np.ones(y3.shape)
y4 = np.array([27, 28, 28, 29, 30])
x4= np.array([74, 74, 75, 75, 75])
x = np.int64(np.concatenate((x1,x2,x3,x4)))
y = np.int64(np.concatenate((y1,y2,y3,y4)))
x = xr.DataArray(x, dims='points')
y = xr.DataArray(y, dims='points')

# read df

ds = df.so.isel(x=x,y=y)
